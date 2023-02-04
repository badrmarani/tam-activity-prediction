import os

DIR = "enzymes_files_fasta/"
OUTPUT_DIR = "enzymes_outputs/"

if not os.path.exists(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)
 
for file in os.listdir(DIR):
    filename = os.fsdecode(file)


    os.system(f"""
        protein-folding-softwares/colabfold_batch/bin/colabfold_batch \
            --amber \
            --templates \
            --num-recycle 3 \
            --use-gpu-relax \
            {DIR} \
            {OUTPUT_DIR} \
    """)