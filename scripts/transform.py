import argparse
import pandas as pd

def main():
    # -------- argparse 设置 --------
    parser = argparse.ArgumentParser(
        description="Transform minimap2 PAF file into clean TSV (first 12 cols)."
    )

    parser.add_argument(
        "--infile", "-i",
        required=True,
        help="Input PAF file (e.g., syn.self.asm5.paf)"
    )

    parser.add_argument(
        "--outfile", "-o",
        required=True,
        help="Output TSV file (e.g., syn_adjust.tsv)"
    )

    args = parser.parse_args()

    infile = args.infile
    outfile = args.outfile

    print(f"[INFO] Reading PAF: {infile}")

    # 读取 PAF
    df = pd.read_csv(infile, sep="\t", header=None, engine="python")

    # 截取前 12 列
    df = df.iloc[:, :12]

    # 统一列名
    df.columns = [
        "qname","qlen","qstart","qend","strand",
        "tname","tlen","tstart","tend","n_match","aln_len","mapq"
    ]

    # 保存
    df.to_csv(outfile, sep="\t", index=False)

    print(f"[INFO] Saved cleaned TSV to: {outfile}")


if __name__ == "__main__":
    main()
