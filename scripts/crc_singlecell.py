from pathlib import Path
import pandas as pd

PROJECT_DIR = Path("/cfs/earth/scratch/suereser/crc_singlecell")

data_file = PROJECT_DIR / "data" / "dataset.csv"

df = pd.read_csv(data_file)

print(df.head())