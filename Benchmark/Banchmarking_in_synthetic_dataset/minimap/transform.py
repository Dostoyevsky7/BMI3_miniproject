import pandas as pd

# PAF 没有 header、每列 tab 分隔、标签字段可变。
df = pd.read_csv("syn.self.asm5.paf", sep="\t", header=None, engine="python")

# 截取前 12 列
df = df.iloc[:, :12]

df.columns = [
    "qname","qlen","qstart","qend","strand",
    "tname","tlen","tstart","tend","n_match","aln_len","mapq"
]

df.to_csv("syn_adjust.tsv", sep="\t", index=False)