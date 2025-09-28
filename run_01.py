# %%
import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt 
import seaborn as sns 
import warnings 
import pathlib 
import os 
import ot
from helper_functions import * 

warnings.filterwarnings("ignore")

from tqdm import tqdm

project = "ecd_wgs_enriched_features"
inputdir = f"/Volumes/HNSD01/storage/{project}"
outputdir = f"/Volumes/HNSD01/outdir/{project}"

path_to_01_output = os.path.join(outputdir, "01_output")
os.system(f"mkdir -p {path_to_01_output}")

feature_name = "FLEN"

all_input_files = [item for item in pathlib.Path(inputdir).glob(f"*/*/*/*/*_GWfeature_*.csv")]
metadata = pd.DataFrame(
    data = all_input_files, columns = ["path"]
)
metadata["SampleID"] = metadata["path"].apply(lambda x: x.name.split("_")[0].split("-")[1])
metadata["Enrich_strategy"] = metadata["path"].apply(lambda x: str(x).split("/")[-5])
metadata["Enrich_panel"] = metadata["path"].apply(lambda x: str(x).split("/")[-4])
metadata["feature_type"] = metadata["path"].apply(lambda x: x.name.split("_")[-1].replace(".csv", ""))
metadata["feature_dir"] = metadata["path"].apply(lambda x: str(x).split("/fragmentomics_features")[0])
mainsrc = "/Users/hieunguyen/src/ecd_wgs_enriched_features"
sample_metadata = pd.read_csv(os.path.join(mainsrc, "metadata.csv"))[["LABCODE", "Label"]]
metadata = metadata.merge(sample_metadata, right_on = "LABCODE", left_on = "SampleID")

for enrich_strategy in metadata.Enrich_strategy.unique().tolist():
    for enrich_panel in metadata[metadata.Enrich_strategy == enrich_strategy].Enrich_panel.unique().tolist():
        print(f"Processing {enrich_strategy} - {enrich_panel}")
        path_to_save_output = os.path.join(path_to_01_output, f"{enrich_strategy}_{enrich_panel}")
        os.system(f"mkdir -p {path_to_save_output}")

        if os.path.isfile(os.path.join(path_to_save_output, "status.txt")) == False:
            run_analysis(
                all_input_samples = metadata[(metadata.Enrich_panel == enrich_panel) & (metadata.Enrich_strategy == enrich_strategy)].feature_dir.unique().tolist(), 
                path_to_save_output = path_to_save_output)
        else:
            print(f"Results exist at {path_to_save_output}, skipping...")