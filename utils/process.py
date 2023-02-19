import os
import warnings
from itertools import zip_longest
from urllib.request import urlopen

import pandas as pd
from tqdm import tqdm


def get_smiles_selfies(
    filename: str = "liste_substrats_accepteurs.txt",
    save: bool = False, output: str = None
) -> pd.DataFrame:

    import requests
    import selfies as sf

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


def monomer_to_dimer(dir):
    files = os.listdir(dir)

    for file in tqdm(files):
        full_path = os.path.join(dir, file)
        with open(full_path, "r") as content:
            lines = content.readlines()

            temp = lines[1:]
            lines[-1] = f"{lines[-1][:-1]}:\n"
            for line in temp:
                lines.append(line)

        open(full_path, "w").writelines(lines)


def rename_fasta_files(dir):
    if not os.path.exists(dir):
        return None

    files = os.listdir(dir)
    for file in tqdm(files):
        full_path = os.path.join(dir, file)

        with open(full_path, "r") as f:
            lines = f.readlines()
            lines[0] = f">{file[:-6]}\n"

        open(
            os.path.join(full_path),
            "w",
        ).writelines(lines)


def get_fasta(output_dir):
    URL = "https://rest.uniprot.org/uniprotkb/{}.fasta"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    with open("list_enzymes.txt") as enzymes:
        for enz in tqdm(enzymes):
            enz = enz[:-1]
            with urlopen(URL.format(enz)) as response:
                filename = f"{output_dir}/{enz}.fasta"
                for line in response:
                    line = line.decode()
                    open(filename, "a").write(line)

                    if line.startswith("Error messages"):
                        warnings.warn(
                            f"ERROR: {enz} is not available in the Uniprot dataset")

            if not os.path.exists(filename):
                warnings.warn(f"ERROR: {enz} has been deleted from UniProtKB")


def predict_enzyme_structure(fasta_dir, output_dir):

    def grouper(iterable, n, fillvalue=None):
        args = [iter(iterable)] * n
        return zip_longest(*args, fillvalue=fillvalue)

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    for f1, f2, f3 in grouper(os.listdir(fasta_dir), 3):
        f1 = os.fsdecode(f1)
        f2 = os.fsdecode(f2)
        f3 = os.fsdecode(f3)

        enz_out1 = output_dir+f1[:6]+"/"
        enz_out2 = output_dir+f2[:6]+"/"
        enz_out3 = output_dir+f3[:6]+"/"
        os.makedirs(enz_out1)
        os.makedirs(enz_out2)
        os.makedirs(enz_out3)

        os.system(f"""
            CUDA_VISIBLE_DEVICES=0 ../protein-folding-softwares/colabfold_batch/bin/colabfold_batch \
                --random-seed 123 \
                --amber \
                --templates \
                --num-models 1 \
                --num-recycle 3 \
                --use-gpu-relax \
                --model-type AlphaFold2-multimer-v2 \
                {fasta_dir+f1} \
                {enz_out1} \
            & \
            CUDA_VISIBLE_DEVICES=1 ../protein-folding-softwares/colabfold_batch/bin/colabfold_batch \
                --random-seed 123 \
                --amber \
                --templates \
                --num-models 1 \
                --num-recycle 3 \
                --use-gpu-relax \
                --model-type AlphaFold2-multimer-v2 \
                {fasta_dir+f2} \
                {enz_out2} \
            & \
            CUDA_VISIBLE_DEVICES=2 ../protein-folding-softwares/colabfold_batch/bin/colabfold_batch \
                --random-seed 123 \
                --amber \
                --templates \
                --num-models 1 \
                --num-recycle 3 \
                --use-gpu-relax \
                --model-type AlphaFold2-multimer-v2 \
                {fasta_dir+f3} \
                {enz_out3} \
        """)


if __name__ == "__main__":
    predict_enzyme_structure(
        fasta_dir=None,
        output_dir=None
    )
