import os
import pickle

import pandas as pd
import torch
from process import get_smiles_selfies
from rdkit import Chem
from torch.utils import data as torch_data
from torchdrug import core, data, models, tasks, utils
from torchdrug.core import Registry
from tqdm import tqdm


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
            # out = {"data": self.data, "act": self.targets}
            # with open(self.out_file, "wb") as f:
            #     pickle.dump(out, f, protocol=pickle.HIGHEST_PROTOCOL)
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
            
            out = {"data": self.data, "act": self.targets}
            with open(self.out_file, "wb") as f:
                pickle.dump(out, f, protocol=pickle.HIGHEST_PROTOCOL)
            

    def get_item(self, index):
        item = {
            "graph1": self.data[index][0],
            "graph2": self.data[index][1],
            "act": self.targets["act"][index]
        }

        if self.transform:
            item = self.transform(item)
        return item

dataset = TAM("dataset/test/", processed_file=True)

model = models.ProteinCNN(input_dim=21,
                          hidden_dims=[1024, 1024],
                          kernel_size=5, padding=2, readout="max")

model2 = models.GIN(input_dim=67,
                    hidden_dims=[256, 256, 256, 256],
                    batch_norm=True, short_cut=True, concat_hidden=True)

task = tasks.InteractionPrediction(
    model, model2=model2, task=dataset.tasks,
    criterion="mse", metric=("mae", "rmse", "spearmanr", "r2"),
    normalization=False, num_mlp_layer=3)

    
optimizer = torch.optim.Adam(task.parameters(), lr=1e-4)
solver = core.Engine(
    task,
    train_set=dataset,
    valid_set=dataset,
    test_set=dataset,
    optimizer=optimizer,
    gpus=[0,1],
    batch_size=32,
    logger="wandb",
    log_interval=10,
)

num_epoch=5000
solver.train(num_epoch=num_epoch)
solver.save(f"checkpoint_{num_epoch}.ckt")
solver.evaluate("train")
