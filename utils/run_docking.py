import os
import gzip
import shutil
import argparse


def run_docking(enz_path: str, sub_path: str, output: str) -> None:
    """
    Here we're using the ADFR Suite.
    * Preparing the receptor (i.e. the transaminase) : Creating the PDBQT file containing
    only the polar hydrogyn charges. The reason I've choose this file extension in particular
    is because it .... (check the youtube video you've seen yeseterday for more reference)
    * preparing the ligand (i.e. the substrat)
     
    """
    if not os.path.exists("PMP.sdf"):
        PMP = "CC1=NC=C(C(=C1O)CN)COP(=O)(O)O"
        os.system("obabel -:'{}' -OPMP.sdf --gen3D".format(PMP)) # convert from smiles to sdf
        os.system("obabel PMP.sdf -OPMP.pdbqt --gen3D") # convert from sdf to pdbqt using obabel
        # os.system("mk_prepare_ligand.py -i PMP.sdf -o PMP.pdbqt") # convert from sdf to pdbqt using adfr
    
    if os.path.exists(output):
        shutil.rmtree(output)
    
    # enz = os.path.join(output, "enzymes")
    os.makedirs(output)
    henz = os.path.join(output, "holoenzymes")
    os.makedirs(henz)
    
    list_enz_path = os.listdir(enz_path)
    list_sub_path = os.listdir(sub_path)
    assert len(list_enz_path) == len(list_sub_path), "ERROR"
    
    for enz, sub in zip(list_enz_path, list_sub_path):
        
        # not efficient... but it's working...
        temp_enz = enz.split(".")
        enz_name = temp_enz[0].split("_")[0]
        temp_enz[0] = enz_name
        temp_enz.insert(1, ".temp.")
        temp_enz = "".join(temp_enz)
        temp_enz_path = os.path.join(output, temp_enz)

        # basic implementation of `grep` with the pattern 'ATOM'
        with (
            open(os.path.join(enz_path, enz), "r") as file,
            open(temp_enz_path, "w") as outfile,
        ):  
            for line in file:
                if "ATOM" in line:
                    outfile.write(line)
            # pass

        # obabel
        os.system("obabel {} -O{} -gen3D".format(temp_enz_path, temp_enz_path))
        
        x = os.path.join(output, enz_name+".pdbqt")
        # os.system("prepare_receptor -r {} -o {}".format(temp_enz_path, x)) # preparing the receptor        
        os.system("obabel {} -O{} -gen3D -xr".format(temp_enz_path, x))


        # gnina
        holoenzyme_name = f"holoenzyme_{enz_name}.pdbqt"
        holoenzyme_path = os.path.join(henz, holoenzyme_name)
        
        os.system("gnina -r {} -l {} --autobox_ligand {} --seed 0 -o {} --exhaustiveness {} --scoring {}".format(
            x, # temp_enz_path,
            "PMP.pdbqt",
            x, # temp_enz_path,
            holoenzyme_path,
            32, # os.cpu_count(),
            "vinardo",
        ))
        
        # apparently, older obabel has trouble with gzipped files...
        # with gzip.open(holoenzyme_path, "rb") as f_in:
        #     with open(os.path.join(henz, holoenzyme_name[:-3]), "wb") as f_out:
        #         f_out.write(f_in.read())
               
        # optional: evalute the docked poses with the rmsd
        # save the output of the command by adding `> log.txt` at the end.
        os.system("obrms --firstonly {} {}".format(temp_enz_path, os.path.join(henz, holoenzyme_name[:-3])))

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--enz-path", type=str)
    parser.add_argument("--sub-path", type=str)
    parser.add_argument("--out-path", type=str)
    args = parser.parse_args()
    
    run_docking(*args.__dict__.values())
