from dgl.dataloading import GraphDataLoader
import os

import dgl
import pandas as pd
import requests
import selfies as sf
import torch
from dgllife.utils import smiles_to_complete_graph
from dgl.data import DGLDataset
from graphein.molecule.config import MoleculeGraphConfig
from graphein.molecule.graphs import construct_graph as consrtuct_molecule
from graphein.protein.config import ProteinGraphConfig
from graphein.protein.graphs import construct_graph as consrtuct_protein
from torch_geometric.data import Data, Dataset
from tqdm import tqdm


def get_smiles_selfies(
        filename: str = "liste_substrats_accepteurs.txt",
    save: bool = False, output: str = None
) -> pd.DataFrame:

    def encoder(mol: str, to_selfies: bool = False):
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{mol}/property/CanonicalSMILES/TXT"
        smiles = requests.get(url).text.rstrip()
        if "NotFound" in smiles:
            return None
        elif to_selfies:
            return sf.encoder(smiles)
        return smiles

    df = pd.DataFrame(columns=["name", "smiles", "selfies"])
    with open(filename, "r") as file:
        lines = file.readlines()
        for x in lines:
            x = x[:-1]
            print(x)
            new_row = pd.Series({
                "name": x,
                "smiles": encoder(x), "selfies": encoder(x, True)
            })

            df = pd.concat(
                [df, new_row.to_frame().T],
                ignore_index=True,
            )
    df.drop_duplicates(subset=["names",], inplace=True)
    if save:
        df.to_csv(
            output,
            sep=";",
            index=False,
        )
    return df

class TADataset(DGLDataset):
    def __init__(
        self,
        raw_dir=None,
        enz_config=None,
        mol_config=None,
        url=None,
        save_dir=None,
        force_reload=False,
        verbose=False,
    ):

        if enz_config is None:
            self.enz_config = ProteinGraphConfig()
        else:
            self.enz_config = ProteinGraphConfig(config=enz_config)

        if mol_config is None:
            self.mol_config = MoleculeGraphConfig()
        else:
            self.mol_config = MoleculeGraphConfig(config=mol_config)

        super(TADataset, self).__init__(
            name="ta_dataset",
            raw_dir=raw_dir,
        )

    def process(self):
        self.enz = []
        self.sub = []
        self.labels = []

        data = pd.read_csv(self.raw_dir, sep=";").reset_index()
        for index, row in tqdm(data.iterrows(), total=data.shape[0]):
            act = torch.tensor(row.activity, dtype=torch.float32)
            self.labels.append(act)

            # constructing enzyme graph
            pdb_path = os.path.join("data/raw/enzymes", row.enzyme, row.enzyme+"_relaxed.pdb")
            enz = consrtuct_protein(config=self.enz_config, pdb_path=pdb_path)
            enz = dgl.from_networkx(enz)  # convert to torch.tensor
            self.enz.append(enz)

            # constructing molecular graphs
            sub = consrtuct_molecule(config=self.mol_config, smiles=row.smiles)
            sub = smiles_to_complete_graph(sub)  # convert to torch.tensor
            self.sub.append(sub)

    def __getitem__(self, idx):
        return (
            self.enz[idx],
            self.sub[idx],
            self.labels[idx],
        )

    def __len__(self):
        return len(self.enz)


# test
data_path = "data/raw/ta_dataset_small.csv"
ds = TADataset(
    raw_dir=data_path,
    enz_config=None,
    mol_config=None,
)


dataloader = GraphDataLoader(ds, batch_size=2, shuffle=False)
for epoch in range(100):
    for a, b, labels in dataloader:
        print(a)
        # print(labels)
        break
    break
