import pandas as pd

DATA_PATH = "ta_dataset.csv"
data = pd.read_csv(DATA_PATH, sep=";")

print(data.head())

