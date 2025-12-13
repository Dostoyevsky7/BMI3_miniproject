import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys
import os

###############################################################################
# Argument parsing
###############################################################################
if len(sys.argv) < 2:
    raise ValueError("Usage: python visualize_for_our_software.py <input_tsv> [--out output_dir]")

infile = None
outdir = "."

# parse parameters
args = sys.argv[1:]
i = 0
while i < len(args):
    if args[i] == "--out":
        outdir = args[i+1]
        i += 2
    else:
        infile = args[i]
        i += 1

if infile is None:
    raise ValueError("Error: You must provide an input TSV file.")

print(f"Reading input file: {infile}")

# create output directory
os.makedirs(outdir, exist_ok=True)

# extract prefix for output filenames
prefix = os.path.splitext(os.path.basename(infile))[0]

###############################################################################
# Load data
###############################################################################
df = pd.read_csv(infile, sep="\t")

# only one chromosome is used, comparing with itself
chrom_length = max(df["end1"].max(), df["end2"].max())

###############################################################################
# 1. Dotplot
###############################################################################
plt.figure(figsize=(7,7))

df["length"] = pd.to_numeric(df["length"], errors="coerce").fillna(1)
df["start1"] = pd.to_numeric(df["start1"], errors="coerce")
df["start2"] = pd.to_numeric(df["start2"], errors="coerce")

sizes = df["length"].astype(float) / 500

color_map = df["strand"].map({
    "+": "#2A9D8F",
    "-": "#E76F51"
})

plt.scatter(
    df["start1"], 
    df["start2"],
    s=sizes,
    c=color_map,
    alpha=0.6,
    edgecolors='none'
)

plt.xlabel("Position on chromosome 1")
plt.ylabel("Position on chromosome 1")
plt.title("Dotplot of Long Repeats (Chr1 Self-alignment)")
plt.tight_layout()

dotplot_path = os.path.join(outdir, f"{prefix}.dotplot.png")
plt.savefig(dotplot_path, dpi=300)
plt.close()

print(f"Dotplot saved to {dotplot_path}")

###############################################################################
# 2. Arc Plot
###############################################################################
df = df.sort_values("length")

plt.figure(figsize=(14, 4))

for _, row in df.iterrows():
    x1 = row["start1"]
    x2 = row["start2"]
    xm = (x1 + x2) / 2
    h = abs(x2 - x1) / 2

    t = np.linspace(0, np.pi, 200)
    x = xm + (x1 - xm) * np.cos(t)
    y = h * np.sin(t)

    alpha = 0.15 + 0.85 * (row["length"] / df["length"].max())
    lw = 0.5 + 3 * (row["length"] / df["length"].max())
    color = "#3A6EA5" if row["strand"] == "+" else "#D1495B"

    plt.plot(x, y, color=color, alpha=alpha, linewidth=lw)

plt.xlim(0, df["end1"].max())
plt.ylim(0)
plt.xlabel("Chromosome 1 position")
plt.title("Enhanced Arc Plot of Long Repeats (Chr1)")
plt.tight_layout()

arcplot_path = os.path.join(outdir, f"{prefix}.arcplot.png")
plt.savefig(arcplot_path, dpi=300)
plt.close()

print(f"Arc plot saved to {arcplot_path}")

###############################################################################
# 3. Repeat Density Plot
###############################################################################
window = 50000

df["mid"] = (df["start1"] + df["end1"]) // 2

bins = range(0, chrom_length + window, window)
hist = df.groupby(pd.cut(df["mid"], bins)).size()

plt.figure(figsize=(14, 3))
plt.bar(range(len(hist)), hist, width=1, color="gray")

plt.xlabel(f"Genomic windows (size = {window} bp)")
plt.ylabel("Repeat count")
plt.title("Repeat Density Along Chromosome 1")
plt.tight_layout()

density_path = os.path.join(outdir, f"{prefix}.density.png")
plt.savefig(density_path, dpi=300)
plt.close()

print(f"Density plot saved to {density_path}")

###############################################################################
# 4. Identity vs position plot
###############################################################################
df["start1"] = pd.to_numeric(df["start1"], errors="coerce")
df["identity"] = pd.to_numeric(df["identity"], errors="coerce")
df["length"] = pd.to_numeric(df["length"], errors="coerce").fillna(1)

sizes = df["length"] / 500

color_map = df["strand"].map({
    "+": "#2A9D8F",
    "-": "#E76F51"
})

plt.figure(figsize=(10, 5))

plt.scatter(
    df["start1"],
    df["identity"],
    s=sizes,
    c=color_map,
    alpha=0.6,
    edgecolors='none'
)

plt.xlabel("Position on chromosome 1", fontsize=12)
plt.ylabel("Identity", fontsize=12)
plt.title("Identity vs. Position (Chr1 Long Repeat Alignments)", fontsize=14)

ax = plt.gca()
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

plt.tight_layout()

identity_path = os.path.join(outdir, f"{prefix}.identity.png")
plt.savefig(identity_path, dpi=300)
plt.close()

print(f"Identity plot saved to {identity_path}")
