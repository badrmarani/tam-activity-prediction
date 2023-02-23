import os
from itertools import combinations, product

import pandas as pd
from tqdm import tqdm

output = "data/raw/test.csv"
enz = pd.read_csv("data/raw/enzyme_infos.csv", sep=";")
sub = pd.read_csv("data/raw/substrat_infos.csv", sep=";")
act = pd.read_csv("data/raw/ta_dataset_v3.csv", sep=";").activity.values

lenz = [list(x) for x in enz.values]
lsub = [list(x) for x in sub.values]

A = product(lenz, repeat=1)
B = combinations(lsub, r=1)

all_columns = list(enz.columns)+list(sub.columns)
all_columns.append("activity")
open(output, "a").writelines(";".join(all_columns)+"\n")

for index, (b, a) in tqdm(enumerate(product(B, A)), total=len(act)):
    ab = [str(x) for x in a[0]] + [str(x) for x in b[0]]
    ab.append(str(act[index]))

    open(output, "a").writelines(
        ";".join(ab) + "\n"
    )

df = pd.read_csv("data/raw/ta_dataset_v4.csv", sep=";")

donors = ["L-ACS", "L-Glutamate", "D-ACS"]

for donor in donors:
    temp = df[df.donor == donor]
    temp.to_csv(
        os.path.join("data/raw/", "ta_dataset_v4_"+donor.replace("-", "_"))+".csv",
        sep=";",
        index=False,
    )
