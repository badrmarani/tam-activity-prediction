import os
import re

import deepchem as dc
import pandas as pd
import rdkit
import requests
import selfies as sf
import torch
from rdkit import Chem
from torch_geometric.data import Data, Dataset
from tqdm import tqdm

DATA_PATH = "data/raw/ta_dataset_v3.csv"
data = pd.read_csv(DATA_PATH, sep=";")


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


class PairData(Data):
    def __inc__(self, key, value, *args):
        if bool(re.search('index_s', key)):
            return self.x_s.size(0)
        if bool(re.search('index_t', key)):
            return self.x_t.size(0)
        else:
            return 0


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
        return [f'ta_dataset_{i}.pt' for i in list(self.data.index)]

    def download(self):
        pass

    def process(self):
        self.data = pd.read_csv(self.raw_paths[0], sep=";").reset_index()
        featurizer = dc.feat.MolGraphConvFeaturizer(use_edges=True)
        for index, row in tqdm(self.data.iterrows(), total=self.data.shape[0]):
            act = torch.tensor(row.activity, dtype=torch.float32)

            # constructing molecular graphs
            mol = Chem.MolFromSmiles(row.smiles)
            f = featurizer._featurize(mol)

            # constructing enzyme graphs
            # not yet implemented

            data = PairData(
                x_s=f.node_features,
                edge_index_s=f.edge_index,
                edge_attr_s=f.edge_features,
                x_t=None,
                edge_index_t=None,
                edge_attr_t=None,
                y=act,
            )

            torch.save(
                data,
                os.path.join(self.processed_dir, f"data_{index}.pt")
            )
            # break

    def _get_label(self, label):
        return torch.tensor([label,], dtype=torch.float32)

    def __len__(self):
        return self.data.size(0)

    def __getitem__(self, index):
        return torch.load(
            os.path.join(self.processed_dir, f"data_{index}.pt")
        )


# test
ds = TADataset(
    root="data/",
    filename="ta_dataset_small.csv",
)
