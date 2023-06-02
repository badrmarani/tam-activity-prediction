from torchdrug import data
from torch.utils import data as torch_data
from torchdrug import data, utils
from torchdrug.core import Registry

from rdkit import Chem
import os

import pandas as pd

from tqdm import tqdm 

from process import get_smiles_selfies
import pickle

@Registry.register("datasets.TAM")
class TAM(data.ProteinDataset):

    ta_dataset = "dataset/ta_dataset_v4_L_ACS.csv"
    out_file = "dataset/xxx_gluta.pkl"
    splits = ["train", "valid", "test_fold", "test_family", "test_superfamily"]
    target_fields = ["activity"]

    def __init__(self, path, processed_file=None, verbose=1, **kwargs):
        path = os.path.expanduser(path)
        if not os.path.exists(path):
            os.makedirs(path)
        self.path = path

        self.transform = None
        self.ta = pd.read_csv(self.ta_dataset, sep=";")
        self.ta = self.ta[(self.ta.induction==1) & (self.ta.expression==1)]
        
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
    
    def split(self):
        keys = ["train", "valid", self.test_split]
        offset = 0
        splits = []
        for split_name, num_sample in zip(self.splits, self.num_samples):
            if split_name in keys:
                split = torch_data.Subset(self, range(offset, offset + num_sample))
                splits.append(split)
            offset += num_sample
        return splits


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
            "enz": self.data[index][0],
            "sub": self.data[index][1],
            "act": self.targets["act"][index]
        }

        if self.transform:
            item = self.transform(item)
        return item
