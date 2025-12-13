import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import argparse
import os

# ------------------------------
# Load and preprocess MUMmer show-coords TSV
# ------------------------------
def load_df(filename):
    cols = [
        "S1", "E1", "S2", "E2",
        "LEN1", "LEN2", "IDY",
        "LENR", "LENQ", "COVR", "COVQ",
        "TAG1", "TAG2"
    ]

    df = pd.read_csv(filename, sep="\t", skiprows=4, names=cols)

    df = df.rename(columns={
        "S1": "start1",
        "E1": "end1",
        "S2": "start2",
        "E2": "end2",
        "IDY": "identity",
        "LEN1": "length"
    })

    df["strand"] = np.where(df["start2"] < df["end2"], "+", "-")
    return df


# ------------------------------
# Argument parsing
# ------------------------------
parser = argparse.ArgumentParser()
parser.add_argument("tsvfile", help="MUMmer show-coords TSV file (e.g., chr1_out.tsv)")
parser.add_argument("--out", default="mummer_plots",
                    help="Directory to store output figures")
args = parser.parse_args()

filename = args.tsvfile
outdir = args.out
os.makedirs(outdir, exist_ok=True)


# ===================================================
# 1. Dotplot
# ===================================================
df = load_df(filename)
chrom_length = max(df["end1"].max(), df["end2"].max())

df["length"] = pd.to_numeric(df["length"], errors="coerce").fillna(1)
df["start1"] = pd.to_numeric(df["start1"], errors="coerce")
df["start2"] = pd.to_numeric(df["start2"], errors="coerce")

plt.figure(figsize=(7, 7))
sizes = df["length"] / 500
color_map = df["strand"].map({"+": "#2A9D8F", "-": "#E76F51"})

plt.scatter(df["start1"], df["start2"], s=sizes, c=color_map,
            alpha=0.6, edgecolors="none")

plt.xlabel("Position on Chromosome 1")
plt.ylabel("Position on Chromosome 1")
plt.title("Dotplot of Chr1 Self-Alignment (MUMmer)")
plt.tight_layout()
plt.savefig(os.path.join(outdir, "dotplot.png"), dpi=300)
plt.close()


# ===================================================
# 2. Arc Plot
# ===================================================
df = load_df(filename)
df = df.sort_values("length")

plt.figure(figsize=(14, 4))
for _, row in df.iterrows():
    x1, x2 = row["start1"], row["start2"]
    xm = (x1 + x2) / 2
    r = abs(x2 - x1) / 2

    t = np.linspace(0, np.pi, 200)
    x = xm + (x1 - xm) * np.cos(t)
    y = r * np.sin(t)

    norm_len = row["length"] / df["length"].max()
    alpha = 0.1 + 0.9 * norm_len
    lw = 0.5 + 3.0 * norm_len
    color = "#3A6EA5" if row["strand"] == "+" else "#D1495B"

    plt.plot(x, y, color=color, alpha=alpha, linewidth=lw)

plt.xlim(0, df["end1"].max())
plt.ylim(0)
plt.xlabel("Chromosome 1 Position")
plt.title("Arc Plot of Long Repeats (Chr1 Self-alignment, MUMmer)")
plt.tight_layout()
plt.savefig(os.path.join(outdir, "arcplot.png"), dpi=300)
plt.close()


# ===================================================
# 3. Density Plot
# ===================================================
df = load_df(filename)
chrom_length = df["end1"].max()
df["mid"] = (df["start1"] + df["end1"]) // 2

window = 50000
bins = np.arange(0, chrom_length + window, window)
hist = df.groupby(pd.cut(df["mid"], bins)).size()

plt.figure(figsize=(14, 3))
plt.bar(range(len(hist)), hist, width=1, color="gray")
plt.xlabel(f"Genomic windows (window size = {window} bp)")
plt.ylabel("Repeat count")
plt.title("Repeat Density Along Chromosome 1 (Self-alignment)")
plt.tight_layout()
plt.savefig(os.path.join(outdir, "density.png"), dpi=300)
plt.close()


# ===================================================
# 4. Identity vs Position
# ===================================================
df = load_df(filename)

df["start1"] = pd.to_numeric(df["start1"], errors="coerce")
df["identity"] = pd.to_numeric(df["identity"], errors="coerce")
df["length"] = pd.to_numeric(df["length"], errors="coerce").fillna(1)

sizes = df["length"] / 500
color_map = df["strand"].map({"+": "#2A9D8F", "-": "#E76F51"})

plt.figure(figsize=(10, 5))
plt.scatter(df["start1"], df["identity"], s=sizes, c=color_map,
            alpha=0.6, edgecolors='none')

plt.ylabel("Identity (%)")
plt.xlabel("Position on chromosome 1")
plt.title("Identity vs. Position (Chr1 Self-alignment)")

ax = plt.gca()
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

plt.tight_layout()
plt.savefig(os.path.join(outdir, "identity.png"), dpi=300)
plt.close()


print(f"All MUMmer plots saved to: {outdir}")