import os
from tqdm import tqdm
from urllib.request import urlopen

URL = "https://rest.uniprot.org/uniprotkb/{}.fasta"
OUTPUT_DIR = "enzymes_fasta/"

if not os.path.exists(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)

with open("list_enzymes_names.txt") as enzymes:
    for enz in tqdm(enzymes, ascii=True):
        enz = enz[:-1]
        with urlopen(URL.format(enz)) as response:
            for line in response:
                line = line.decode()
                open(f"{OUTPUT_DIR}/{enz}.fasta", "a").write(line)
        # break