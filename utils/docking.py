import os

enz_name = "C7RAM1"
sub_name = "2OATA"
sub_smiles_seq = "C(=O)C(=O)O"

# Convert SMILES to SDF
if not os.path.exists(sub_name+".sdf"):
    open(sub_name+".smi", "w").writelines(sub_smiles_seq+"\n")

    os.system(f"""
        obabel \
        {sub_name}.smi \
        -O {sub_name}.sdf \
        --gen3d
    """)

# Docking with GNINA
os.system(f"""
    gnina/build/bin/gnina \
    -r {enz_name}.pdb \
    -l {sub_name}.sdf \
    --autobox_ligand {enz_name}.pdb \
    -o docked_{enz_name.lower()}_{sub_name.lower()}.sdf.gz
    --exhaustiveness 64
""")
