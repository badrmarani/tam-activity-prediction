import parmed as pmd
import argparse


def run(pdb1, pdb2, outp):
    pdb1 = pmd.load_file(pdb1)
    pdb2 = pmd.load_file(pdb2)

    combined = pdb1 + pdb2
    combined.save(outp)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--pdb1", type=str)
    parser.add_argument("--pdb2", type=str)
    parser.add_argument("--outp", type=str)
    args = parser.parse_args()
    
    run(*args.__dict__.values())