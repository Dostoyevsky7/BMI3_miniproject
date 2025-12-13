import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse
import os

# ------------------------------
# 0. Load PAF and save TSV to SAME DIR as original PAF
# ------------------------------
def load_paf(path, save_tsv=True):
    df = pd.read_csv(path, sep="\t", header=None, engine="python")

    df = df.iloc[:, :12]
    df.columns = [
        "qname", "qlen", "qstart", "qend", "strand",
        "tname", "tlen", "tstart", "tend",
        "n_match", "aln_len", "mapq"
    ]

    df = df.rename(columns={
        "qstart": "start1",
        "qend":   "end1",
        "tstart": "start2",
        "tend":   "end2",
        "aln_len": "length"
    })

    numeric_cols = ["start1","end1","start2","end2",
                    "length","qlen","tlen","n_match"]
    for col in numeric_cols:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    chrom_length = int(max(df["qlen"].max(), df["tlen"].max()))
    df["identity"] = df["n_match"] / df["length"] * 100

    # TSV 保存路径：与原 PAF 同目录，且保持原有命名逻辑
    if save_tsv:
        paf_basename = os.path.splitext(path)[0]        # e.g. results/arabidopsis/chr1.self.minimap
        tsv_path = paf_basename + "_transform.tsv"       # → results/arabidopsis/chr1.self.minimap_transform.tsv
        df.to_csv(tsv_path, sep="\t", index=False)
        print(f"PAF processed and saved to: {tsv_path}")

    return df, chrom_length


# ------------------------------
# 1. Dotplot
# ------------------------------
def plot_dotplot(df, chrom_length, out_file):
    plt.figure(figsize=(7, 7))

    sizes = df["length"] / 500
    colors = df["strand"].map({"+": "#2A9D8F", "-": "#E76F51"})

    plt.scatter(df["start1"], df["start2"], s=sizes, c=colors,
                alpha=0.6, edgecolors="none")

    plt.xlim(0, chrom_length)
    plt.ylim(0, chrom_length)
    plt.xlabel("Query position")
    plt.ylabel("Target position")
    plt.title("Dotplot (minimap2 self-alignment)")
    plt.tight_layout()
    plt.savefig(out_file, dpi=300)
    plt.close()


# ------------------------------
# 2. Arc plot
# ------------------------------
def plot_arc(df, chrom_length, out_file, min_len=2000):
    df2 = df[df["length"] >= min_len]
    df2 = df2[(df2["start1"] - df2["start2"]).abs() > 1000]
    df2 = df2.sort_values("length")

    plt.figure(figsize=(14, 4))

    for _, r in df2.iterrows():
        x1, x2 = r["start1"], r["start2"]
        xm = (x1 + x2) / 2
        r0 = abs(x2 - x1) / 2

        t = np.linspace(0, np.pi, 150)
        x = xm + (x1 - xm) * np.cos(t)
        y = r0 * np.sin(t)

        color = "#3A6EA5" if r["strand"] == "+" else "#D1495B"
        alpha = 0.1 + 0.9 * r["length"] / df2["length"].max()
        lw = 0.4 + 2.5 * r["length"] / df2["length"].max()

        plt.plot(x, y, color=color, alpha=alpha, linewidth=lw)

    plt.xlim(0, chrom_length)
    plt.ylim(0)
    plt.xlabel("Chromosome position")
    plt.title("Arc plot (minimap2 repeats)")
    plt.tight_layout()
    plt.savefig(out_file, dpi=300)
    plt.close()


# ------------------------------
# 3. Density
# ------------------------------
def plot_density(df, chrom_length, out_file, window=50000):
    df["mid"] = (df["start1"] + df["end1"]) // 2
    bins = np.arange(0, chrom_length + window, window)
    hist = df.groupby(pd.cut(df["mid"], bins)).size()

    plt.figure(figsize=(14, 3))
    plt.bar(range(len(hist)), hist, width=1, color="gray")
    plt.xlabel(f"Genome windows ({window} bp)")
    plt.ylabel("Repeat count")
    plt.title("Repeat density along chromosome")
    plt.tight_layout()
    plt.savefig(out_file, dpi=300)
    plt.close()


# ------------------------------
# 4. Identity plot
# ------------------------------
def plot_identity(df, out_file):
    sizes = df["length"] / 500
    colors = df["strand"].map({"+": "#2A9D8F", "-": "#E76F51"})

    plt.figure(figsize=(10, 4))
    plt.scatter(df["start1"], df["identity"], s=sizes,
                c=colors, alpha=0.6, edgecolors="none")

    plt.xlabel("Position on genome")
    plt.ylabel("Identity (%)")
    plt.title("Identity vs Position")
    plt.tight_layout()
    plt.savefig(out_file, dpi=300)
    plt.close()


# ------------------------------
#  MAIN
# ------------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("paf_file", help="PAF file from minimap output")
    parser.add_argument("--out", default="minimap_plots",
                        help="Directory to store output figures")
    args = parser.parse_args()

    input_file = args.paf_file
    outdir = args.out
    os.makedirs(outdir, exist_ok=True)

    df, chrom_len = load_paf(input_file, save_tsv=True)

    plot_dotplot(df, chrom_len, os.path.join(outdir, "dotplot.png"))
    plot_arc(df, chrom_len, os.path.join(outdir, "arcplot.png"))
    plot_density(df, chrom_len, os.path.join(outdir, "density.png"))
    plot_identity(df, os.path.join(outdir, "identity.png"))

    print(f"All plots saved under: {outdir}")
