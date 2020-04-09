import pandas as pd


# it first requires to remove the string 'PDB' before the IDs in the files
# that is the only position with 1 sigle space and it is not parsed correctly
# can be done easier in the shell like this:
# sed -i 's/PDB / /g' ./pdbfold_*.dat


def get_pdbfold_df(filepath):
    with open(filepath) as dat_filein:
        df = pd.read_csv(dat_filein, skiprows=(0, 1, 2, 3), sep="\s+")
    for i in range(
        len(df["Target"])
    ):  # this slpits the ID:chain in target for easier processing further
        df["Target"][i] = df["Target"][i].split(":")
    return df
