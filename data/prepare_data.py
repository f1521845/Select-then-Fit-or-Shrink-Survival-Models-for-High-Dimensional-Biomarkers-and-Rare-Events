"""Export the SurvSet datasets used in the real-data notebooks to CSV.

SurvSet stores its data as Python pickles, so the R notebooks call this script
once (it is skipped when the CSVs already exist):

    pip install SurvSet
    python data/prepare_data.py

Writes data/support2.csv, data/DLBCL.csv and data/gse1992.csv. PBC comes from
the R survival package and CHOP/R-CHOP from the bujar package data, so they are
loaded directly in their notebooks.
"""
from pathlib import Path

from SurvSet.data import SurvLoader

OUT = Path(__file__).resolve().parent
DATASETS = ["support2", "DLBCL", "gse1992"]


def main() -> None:
    loader = SurvLoader()
    for name in DATASETS:
        path = OUT / f"{name}.csv"
        if path.exists():
            print(f"{path.name}: exists, skipped")
            continue
        df, ref = loader.load_dataset(name).values()
        df.to_csv(path, index=False)
        print(f"{path.name}: {df.shape[0]} rows, {df.shape[1] - 3} predictors, "
              f"{int(df['event'].sum())} events  (source: {ref})")


if __name__ == "__main__":
    main()
