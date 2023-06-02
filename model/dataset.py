import os
import pickle

import dgl
import pandas as pd
import torch
from biopandas.pdb import PandasPdb
from graphein.ml.conversion import GraphFormatConvertor
from graphein.protein.config import ProteinGraphConfig
from graphein.protein.graphs import construct_graph
from process import get_smiles_selfies
from rdkit import Chem
from torch.utils import data as torch_data
from torch_geometric.data import Data, Dataset
from torch_geometric.utils import smiles
from torchdrug import data, utils
from torchdrug.core import Registry
from tqdm import tqdm


def pdb_to_graph(filename, contain_b_factor=True):

    def get_distance(coords):
        diff = coords.unsqueeze(1) - coords.unsqueeze(0)
        distance = torch.sqrt(diff.pow(2).sum(-1))
        return distance

    atom_df = PandasPdb().read_pdb(filename)
    atom_df = atom_df.df["ATOM"]

    aa_df = atom_df.groupby("residue_number", as_index=False)[
        ["x_coord", "y_coord", "z_coord", "b_factor"]
    ].mean().sort_values("residue_number")

    coords = torch.from_numpy(aa_df[["x_coord", "y_coord", "z_coord"]].values)
    distance = get_distance(coords)
    adj = distance < 6.0
    u = torch.nonzero(adj)[:, 0]
    v = torch.nonzero(adj)[:, 1]
    graph = dgl.graph((u, v), num_nodes=len(coords))

    if contain_b_factor:
        b_factor = torch.from_numpy(aa_df["b_factor"].values)
        graph.ndata["b_factor"] = b_factor
    return graph


@Registry.register("datasets.TAM")
class TAM(data.ProteinDataset):

    ta_dataset = "dataset/ta_dataset_v4_L_Glutamate.csv"
    out_file = "dataset/tttam_l_glutamate.pkl"
    target_fields = ["activity"]

    def __init__(self, path, processed_file=None, verbose=1, **kwargs):
        path = os.path.expanduser(path)
        if not os.path.exists(path):
            os.makedirs(path)
        self.path = path

        self.transform = None
        self.ta = pd.read_csv(self.ta_dataset, sep=";")
        self.data = []
        self.targets = {"act": []}

        if not processed_file:
            self.load_pairs(verbose=verbose, **kwargs)
            out = {"data": self.data, "act": self.targets}
            with open(self.out_file, "wb") as f:
                pickle.dump(out, f, protocol=pickle.HIGHEST_PROTOCOL)
        else:
            with open(self.out_file, "rb") as f:
                out = pickle.load(f)
            self.data = out["data"]
            self.targets = out["act"]
    

    def load_pairs(self, transform=None, lazy=None, verbose=0):
        self.transform = transform
        self.lazy = lazy

        df = tqdm(self.ta.iterrows(), total=self.ta.shape[0]) if verbose else self.ta.iterrows()
        for idx, row in df:
            pdb_file = os.path.join("dataset/enzymes", row.enzyme, f"{row.enzyme}_relaxed.pdb")
            graph1 = data.Protein.from_pdb(
                pdb_file=pdb_file,
                atom_feature="default",
                bond_feature="default",
                residue_feature="symbol",
                mol_feature=None,
                kekulize=True)
            
            seq = get_smiles_selfies(row.acceptor_name)["smiles"]
            mol = Chem.MolFromSmiles(seq)
            graph2 = data.Molecule.from_molecule(
                mol,
                atom_feature="default",
                bond_feature="default",
                kekulize=True)

            self.data.append([graph1, graph2])
            self.targets["act"].append(row.activity)
            
            # out = {"data": self.data, "act": self.targets}
            # with open(self.out_file, "wb") as f:
            #     pickle.dump(out, f, protocol=pickle.HIGHEST_PROTOCOL)
            

    def get_item(self, index):
        item = {
            "graph1": self.data[index][0],
            "graph2": self.data[index][1],
            "act": self.targets["act"][index]
        }

        if self.transform:
            item = self.transform(item)
        return item


class PairData(Data):
    def __init__(self,
        x_a=None, edge_index_a=None, edge_attr_a=None,
        x_b=None, edge_index_b=None, edge_attr_b=None,
        y=None, pos=None, **kwargs,
    ) -> None:
        super().__init__()

        # enzyme graph
        self.x_a = x_a
        self.edge_index_a = edge_index_a
        self.edge_attr_a = edge_attr_a

        # substrat graph
        self.x_b = x_b
        self.edge_index_b = edge_index_b
        self.edge_attr_b = edge_attr_b
        self.y = y

    def __inc__(self, key, value, *args, **kwargs):
        if key == "edge_index_a":
            return self.x_a.size(0)
        if key == "edge_index_b":
            return self.x_b.size(0)
        else:
            return super().__inc__(key, value, *args, **kwargs)


class TADataset(Dataset):
    def __init__(self, root, filename, transform=None, pre_transform=None) -> None:
        self.filename = filename
        super(TADataset, self).__init__(root, transform, pre_transform)

    @property
    def raw_file_names(self):
        return self.filename

    @property
    def processed_file_names(self):
        self.data = pd.read_csv(self.raw_paths[0], sep=";").reset_index()
        return [f"ta_dataset_{i}.pt" for i in list(self.data.index)]

    def download(self):
        pass

    def process(self):
        self.all_data = []
        self.data = pd.read_csv(self.raw_paths[0], sep=";").reset_index()
        for index, row in tqdm(self.data.iterrows(), total=self.data.shape[0]):
            act = torch.tensor(row.activity, dtype=torch.float32)

            # constructing molecular graphs
            ga = smiles.from_smiles(row.smiles)

            # constructing enzyme graphs
            config = ProteinGraphConfig()
            pdb_path = os.path.join(self.raw_paths[1], row.enzyme+"_relaxed.pdb")
            gb = construct_graph(config=config, pdb_path=pdb_path)
            gb = GraphFormatConvertor(
                "nx", "pyg",
                verbose="gnn",
                columns=None
            ).convert_nx_to_pyg(g)

            data = PairData(
                ga.x, ga.edge_index, ga.edge_attr,  # substrat
                None, gb.edge_index, None,          # enzyme
                y=act,
            )
            self.all_data.append(data)
        torch.save(
            self.all_data,
            os.path.join(self.processed_dir, f"data.pt")
        )

    def _get_label(self, label):
        return torch.tensor([label,], dtype=torch.float32)

    def __len__(self):
        return self.data.shape[0]

    def __getitem__(self, index):
        data = torch.load(
            os.path.join(self.processed_dir, f"data.pt")
        )
        return data[index]
