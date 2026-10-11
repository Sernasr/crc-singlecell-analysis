"""Explore CRC single-cell datasets and inspect their metadata and expression matrices."""

from pathlib import Path
from io import StringIO

import matplotlib
matplotlib.use("Agg")  # Save figures on a remote server without a display
import matplotlib.pyplot as plt
import pandas as pd
from bs4 import BeautifulSoup

# 1. Setup

PROJECT_DIR = Path("/cfs/earth/scratch/suereser/crc_singlecell")
DATA_DIR = PROJECT_DIR / "data"
RESULTS_DIR = PROJECT_DIR / "results" / "exploration"
GSE108989_RESULTS = RESULTS_DIR / "GSE108989"
GSE146771_RESULTS = RESULTS_DIR / "GSE146771"

for folder in (GSE108989_RESULTS, GSE146771_RESULTS):
    folder.mkdir(parents=True, exist_ok=True)

GSE108989_DIR = DATA_DIR / "2018_GSE108989"
GSE146771_DIR = DATA_DIR / "2020_GSE146771"

# 2. Helper Functions

def list_files(folder):
    """
    List all files in a folder together with their size in MB.
    """
    for file in sorted(folder.iterdir()):
        if file.is_file():
            size_mb = file.stat().st_size / (1024 ** 2)
            print(f"{file.name} | {size_mb:.2f} MB")

def inspect_table(path, sep="\t", nrows=5):
    """
    Load only the first rows of a compressed or uncompressed table
    to inspect its structure without loading the full dataset.
    """
    df = pd.read_csv(
        path,
        sep=sep,
        nrows=nrows,
        compression="infer"
    )

    print(f"File: {path.name}")
    print(f"Preview shape: {df.shape}")
    print(f"Number of columns: {len(df.columns)}")
    print("\nColumns:")
    print(df.columns.tolist())

    print("\nFirst rows:")
    print(df.head())

    return df

def inspect_whitespace_table(path, nrows=5):
    """
    Inspect tables that are separated by whitespace rather than tabs.
    """
    df = pd.read_csv(
        path,
        sep=r"\s+",
        nrows=nrows,
        compression="infer"
    )

    print(f"File: {path.name}")
    print(f"Preview shape: {df.shape}")
    print(f"Number of columns: {len(df.columns)}")

    print("\nFirst column names:")
    print(df.columns[:10].tolist())

    print("\nFirst rows:")
    print(df.iloc[:, :10])

    return df

# 3. GSE108989, Zhang 2018

# 3.1 Available files

list_files(GSE108989_DIR)

# 3.2 Inspect GSE108989 files

gse108989_mutation = inspect_table(
    GSE108989_DIR / "GSE108989_CRC.bulk.exome.mutation.tab.gz"
)

gse108989_bulk_count = inspect_table(
    GSE108989_DIR / "GSE108989_CRC.bulk.S12.count.tab.gz"
)

gse108989_bulk_tpm = inspect_table(
    GSE108989_DIR / "GSE108989_CRC.bulk.S12.TPM.tab.gz"
)

# 3.3 Cell metadata

html_path = GSE108989_DIR / "GSE108989_cell_metadata_full.htm"

with open(html_path, "r", encoding="utf-8") as f:
    soup = BeautifulSoup(f, "html.parser")

text = soup.find("pre").get_text()

gse108989_metadata = pd.read_csv(
    StringIO(text),
    sep="\t",
    comment="#"
)

print("Shape:", gse108989_metadata.shape)
print("Columns:", gse108989_metadata.columns.tolist())

print(gse108989_metadata.head())
gse108989_metadata.columns = gse108989_metadata.columns.str.strip()
print(gse108989_metadata.columns.tolist())

gse108989_count_preview = pd.read_csv(
    GSE108989_DIR / "GSE108989_CRC.TCell.S11138.count.txt.gz",
    sep="\t",
    nrows=5,
    compression="infer"
)

# 3.3.1 Metadata overview

print("Number of cells:", gse108989_metadata.shape[0])
print("Number of patients:", gse108989_metadata["Patient_ID"].nunique())
print("Number of clusters:", gse108989_metadata["majorCluster"].nunique())
print("Number of sample types:", gse108989_metadata["sampleType"].nunique())

print("Cells per patient:")
print(gse108989_metadata["Patient_ID"].value_counts())

print("Major clusters:")
print(gse108989_metadata["majorCluster"].value_counts())

print("Sample types:")
print(gse108989_metadata["sampleType"].value_counts())

count_cells = set(gse108989_count_preview.columns[2:])
metadata_cells = set(gse108989_metadata["UniqueCell_ID"])

print("Cells in count preview:", len(count_cells))
print("Metadata cells:", len(metadata_cells))
print("Matching cells:", len(count_cells & metadata_cells))

# 3.4 Single cell expression matrices

gse108989_count_preview = inspect_table(
    GSE108989_DIR / "GSE108989_CRC.TCell.S11138.count.txt.gz",
    nrows=5
)

gse108989_tpm_preview = inspect_table(
    GSE108989_DIR / "GSE108989_CRC.TCell.S11138.TPM.txt.gz",
    nrows=5
)

gse108989_norm_preview = inspect_table(
    GSE108989_DIR / "GSE108989_CRC.TCell.S10805.norm.centered.txt.gz",
    nrows=5
)

# 4. GSE146771, 2020

# 4.1 Available files

list_files(GSE146771_DIR)

# 4.2 Smart-seq2 metadata

smartseq_metadata = inspect_table(
    GSE146771_DIR / "GSE146771_CRC.Leukocyte.Smart-seq2.Metadata.txt.gz"
)

smartseq_meta = pd.read_csv(
    GSE146771_DIR / "GSE146771_CRC.Leukocyte.Smart-seq2.Metadata.txt.gz",
    sep="\t"
)

print("Number of cells:", smartseq_meta.shape[0])
print("Number of samples:", smartseq_meta["Sample"].nunique())
print("Number of global cell types:", smartseq_meta["Global_Cluster"].nunique())
print("Number of subclusters:", smartseq_meta["Sub_Cluster"].nunique())

print("Global cell types:")
print(smartseq_meta["Global_Cluster"].value_counts())

print("Tissue types:")
print(smartseq_meta["Tissue"].value_counts())

# 4.2.1 Dataset overview

print("Number of cells:", smartseq_meta.shape[0])
print("Number of metadata columns:", smartseq_meta.shape[1])
print("Number of samples:", smartseq_meta["Sample"].nunique())

print("Tissue distribution:")
print(smartseq_meta["Tissue"].value_counts())

print("Global cell types:")
print(smartseq_meta["Global_Cluster"].value_counts())

print("Number of subclusters:", smartseq_meta["Sub_Cluster"].nunique())

print(
    smartseq_meta["Sub_Cluster"]
    .value_counts()
    .head(20)
)

print("Cells per sample:")
print(smartseq_meta["Sample"].value_counts())

print("Cell types by tissue:")
print(pd.crosstab(
    smartseq_meta["Global_Cluster"],
    smartseq_meta["Tissue"]
))

# 4.3 10x metadata

tenx_metadata = inspect_table(
    GSE146771_DIR / "GSE146771_CRC.Leukocyte.10x.Metadata.txt.gz"
)

# 4.4 Smart-seq2 TPM matrix

smartseq_tpm_preview = inspect_whitespace_table(
    GSE146771_DIR / "GSE146771_CRC.Leukocyte.Smart-seq2.TPM.txt.gz",
    nrows=5
)

# 4.5 10x TPM matrix

tenx_tpm_preview = inspect_whitespace_table(
    GSE146771_DIR / "GSE146771_CRC.Leukocyte.10x.TPM.txt.gz",
    nrows=5
)


# 5. Data visualization
# All figures are saved in results/exploration/<dataset>/ as PNG files.

def save_barplot(series, output_path, title, xlabel, ylabel,
                 horizontal=False, figsize=(10, 6)):
    """Plot category counts and save the figure."""
    fig, ax = plt.subplots(figsize=figsize)
    series.plot(kind="barh" if horizontal else "bar", ax=ax)
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if not horizontal:
        ax.tick_params(axis="x", labelrotation=45)
        for label in ax.get_xticklabels():
            label.set_ha("right")
    fig.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved figure: {output_path}")


# 5.1 GSE108989 - Cells per patient
save_barplot(
    gse108989_metadata["Patient_ID"].value_counts(),
    GSE108989_RESULTS / "cells_per_patient.png",
    "GSE108989 - Cells per Patient",
    "Patient ID", "Number of Cells"
)

# 5.2 GSE108989 - Cell clusters
save_barplot(
    gse108989_metadata["majorCluster"].value_counts().sort_values(),
    GSE108989_RESULTS / "cell_cluster_distribution.png",
    "GSE108989 - Major Cluster Distribution",
    "Number of Cells", "Major Cluster", horizontal=True
)

# 5.3 GSE108989 - Sample types
save_barplot(
    gse108989_metadata["sampleType"].value_counts(),
    GSE108989_RESULTS / "sample_type_distribution.png",
    "GSE108989 - Sample Types",
    "Sample Type", "Number of Cells"
)

# 5.4 GSE146771 Smart-seq2 - Cells per sample
# Metadata includes Sample, but a separate Patient_ID column is not assumed.
save_barplot(
    smartseq_meta["Sample"].value_counts(),
    GSE146771_RESULTS / "cells_per_sample.png",
    "GSE146771 Smart-seq2 - Cells per Sample",
    "Sample", "Number of Cells", figsize=(12, 6)
)

# 5.5 GSE146771 Smart-seq2 - Global cell types
save_barplot(
    smartseq_meta["Global_Cluster"].value_counts().sort_values(),
    GSE146771_RESULTS / "cell_type_distribution.png",
    "GSE146771 Smart-seq2 - Global Cell Types",
    "Number of Cells", "Global Cell Type", horizontal=True
)

# 5.6 GSE146771 Smart-seq2 - Tissue distribution
save_barplot(
    smartseq_meta["Tissue"].value_counts(),
    GSE146771_RESULTS / "tissue_distribution.png",
    "GSE146771 Smart-seq2 - Tissue Distribution",
    "Tissue", "Number of Cells"
)

print("\nDataset exploration complete. Figures saved under:", RESULTS_DIR)



# ============================================================
# MSI/MSS ANALYSIS AND RESULTS
# ============================================================

# Published MSI-H patient IDs
msi_2018 = {"P0413", "P0825", "P0123", "P0909"}
msi_2020 = {"P0123", "P0413", "P0825"}

# Published MSS patient IDs for GSE108989
mss_2018 = {
    "P0215", "P0411", "P0701", "P1012",
    "P1207", "P1212", "P1228", "P0309"
}

# Verify the 10 Smart-seq2 patient IDs
expected_2020 = {
    "P0104", "P0305", "P0309", "P0411", "P0413",
    "P0720", "P0728", "P0825", "P1212", "P1228"
}

# Load cell metadata for MSI/MSS analysis
from bs4 import BeautifulSoup
from io import StringIO

with open(
    GSE108989_DIR / "GSE108989_cell_metadata_full.htm",
    encoding="utf-8"
) as f:
    soup = BeautifulSoup(f, "html.parser")

meta_2018 = pd.read_csv(
    StringIO(soup.find("pre").get_text()),
    sep="\t",
    comment="#"
)
meta_2018.columns = meta_2018.columns.str.strip()

meta_2020 = pd.read_csv(
    GSE146771_DIR / "GSE146771_CRC.Leukocyte.Smart-seq2.Metadata.txt.gz",
    sep="\t"
)

# Validate patient IDs before assigning clinical labels
ids_2018 = set(meta_2018["Patient_ID"].astype(str).str.strip())
ids_2020 = set(meta_2020["Sample"].astype(str).str.strip())

assert ids_2018 == (msi_2018 | mss_2018)
assert ids_2020 == expected_2020

# Assign MSI/MSS status
map_2018 = {
    **{p: "MSI-H" for p in msi_2018},
    **{p: "MSS" for p in mss_2018}
}

meta_2018["MSI_Status"] = meta_2018["Patient_ID"].map(map_2018)
meta_2020["MSI_Status"] = meta_2020["Sample"].apply(
    lambda p: "MSI-H" if p in msi_2020 else "MSS"
)

# Create output directories
results_dir = PROJECT_DIR / "results" / "exploration"
out_2018 = results_dir / "GSE108989"
out_2020 = results_dir / "GSE146771"

out_2018.mkdir(parents=True, exist_ok=True)
out_2020.mkdir(parents=True, exist_ok=True)

# Generate group summaries
rows = []

for name, df, patient_col, folder in [
    ("GSE108989", meta_2018, "Patient_ID", out_2018),
    ("GSE146771", meta_2020, "Sample", out_2020)
]:
    patient_labels = df[[patient_col, "MSI_Status"]].drop_duplicates()
    patient_labels.to_csv(
        folder / "msi_patient_annotations.csv",
        index=False
    )

    for status in ["MSI-H", "MSS"]:
        subset = df[df["MSI_Status"] == status]

        rows.append({
            "Dataset": name,
            "MSI_Status": status,
            "Patient_Count": subset[patient_col].nunique(),
            "Cell_Count": len(subset),
            "Cell_Percentage": round(len(subset) / len(df) * 100, 2)
        })

comparison = pd.DataFrame(rows)

# Save group summaries for each dataset
for name, folder in [
    ("GSE108989", out_2018),
    ("GSE146771", out_2020)
]:
    comparison[comparison["Dataset"] == name].to_csv(
        folder / "msi_group_summary.csv",
        index=False
    )

# Save final comparison
comparison.to_csv(
    results_dir / "msi_mss_final_comparison.csv",
    index=False
)

# Plot patient and cell counts
for metric, filename in [
    ("Patient_Count", "msi_mss_patient_comparison.png"),
    ("Cell_Count", "msi_mss_cell_comparison.png")
]:
    comparison.pivot(
        index="Dataset",
        columns="MSI_Status",
        values=metric
    ).plot(kind="bar", figsize=(9, 5), rot=0)

    plt.title(metric.replace("_", " ") + " by MSI/MSS Status")
    plt.ylabel(metric.replace("_", " "))
    plt.xlabel("Dataset")
    plt.tight_layout()
    plt.savefig(results_dir / filename, dpi=300, bbox_inches="tight")
    plt.close()

print("MSI/MSS analysis completed and saved.")