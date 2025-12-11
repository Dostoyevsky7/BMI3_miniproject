# BMI3_miniproject

FindLongRepeat is a code developmed by our group to find long repeat in a certain genome.

The whole pipeline is as follows:

# Environment Setup

This repository contains scripts and tools for visualizing and analyzing genome alignments, 
including Minimap2 PAF files and MUMmer4 show-coords outputs.

To ensure full reproducibility, this document provides the exact development environment 
(WSL2 + Ubuntu + Minimap2 + MUMmer4 + Conda), along with step-by-step instructions for rebuilding it from scratch.

## 1. System Environment (WSL2)

The project was developed on **Windows 11** using **WSL2** running Ubuntu 24.04.

System details:

```bash
Distributor ID: Ubuntu
Description: Ubuntu 24.04.3 LTS
Release: 24.04
Codename: noble
Linux DESKTOP-OEGU2PR 6.6.87.2-microsoft-standard-WSL2 #1 SMP PREEMPT_DYNAMIC Thu Jun 5 18:30:46 UTC 2025 x86_64
```

If you are also using Windows, install WSL2:

```bash
wsl --install
```

## 2. Required External Tools

### Minimap2
```bash
Minimap2 --version
2.26-r1175
```

### MUMmer4
```bash
nucmer --version
4.0.1

mummer --version
4.0.1
```


If these are missing, install via APT:
```bash
sudo apt update
sudo apt install Minimap2 mummer
```

## 3. Conda Environment Setup

All Python dependencies are stored in:
```bash
environment.yml
```

### Create the same environment:
```bash
conda env create -f environment.yml
```

### Activate it:

```bash
conda activate base
```

## 4. System Package List (APT)

This project also includes a full list of installed APT packages:
```bash
apt_installed.txt
```

This file ensures bit-level reproducibility.

To restore all APT packages (optional):

```bash
sudo xargs -a apt_installed.txt apt install -y
```

---

# Repository Contents

| File | Description |
|------|-------------|
| `environment.yml` | Conda environment for all Python dependencies |
| `apt_installed.txt` | Complete list of APT-installed system packages |
| `scripts/` | Analysis and visualization scripts |
| `data/` | Example alignment outputs (optional) |

---

# Pipeline

## 1. generate synthesized data

First, we used data.py to generate a small test dataset by running:
```bash
python scripts/data.py --out data/synthetic/
```
This script creates a synthetic genome containing 20 implanted repeats (synthetic_genome.fasta) along with the corresponding ground-truth list (ground_truth_repeats.tsv).

## 2. find long repeat in the synthesized data

Based on the synthesized data, the core code `repeatfinder.py` could be executed to find long repeat:
```bash
python scripts/repeatfinder.py find \
    --genome data/synthetic/synthetic_genome.fasta \
    --out results/synthetic/far-fast.tsv \
    --chromosome 1 \
    --mask 11110111101111 \
    --step 4 \
    --nseeds 15 \
    --max-occ 300 \
    --band 400 \
    --minlen 20000
```
This execution would generate a .tsv file that contains the long repeat results.
This script identifies repeat regions within the input genome. Given a dataset as input, it produces a TSV file containing the coordinates, lengths and number of repeats, together with a homology score for each repeat.
When running the tool, the user must specify the input file and a set of parameters. These parameters govern the speed, sensitivity, and memory usage of the algorithm. To make the tool more convenient, we provide several predefined parameter presets—such as fast and balanced. Users may either rely on these presets or manually set their own parameter combinations.

## 3. accuracy

`evaluate_repeats.py` is used to evaluate the results of the generated long repeats, comparing to the ground truth synthesized.
```bash
python scripts/evaluate_repeats.py \
    --truth data/synthetic/ground_truth_repeats.tsv \
    --pred results/synthetic/far-fast.tsv \
    --overlap 0.8
```

## 4. Benchmark: other package performance

**Minimap2** and **MUMmer4** are 2 tools rely on **seed-based approach**:

- Minimap2 is a fast, flexible long-read and genome aligner (https://github.com/lh3/minimap2)
- MUMmer4 is a high-precision genome comparison toolkit (https://github.com/mummer4/mummer)

When aligning a genome to itself, any repeated segment will appear as an off-diagonal match, making both tools powerful for detecting long genomic repeats.

Since our algorithm is also **seed-based**, we compare our algorthm with Minimap2 and MUMmer4.

### Minimap2
#### run Minimap2 on the synthesized data:
```bash
minimap2 -x asm5 -t 1 \
-N 1000 \
-p 0.0   \
data/synthetic/synthetic_genome.fasta data/synthetic/synthetic_genome.fasta > results/synthetic/syn.minimap.paf
```
#### Minimap2 accuracy
**Minimap2** generates a `.paf` file as the result. For the convenience to calculate the accuracy, following command could transform the result into a `.tsv` file:
```bash
python scripts/transform.py \
    -i results/synthetic/syn.minimap.paf \
    -o results/synthetic/syn_transform_minimap.tsv
```
After transformed into `.tsv`, the accuracy could be calculated by:
```bash
python scripts/evaluate_repeats.py \
    --truth data/synthetic/ground_truth_repeats.tsv \
    --pred results/synthetic/syn_transform_minimap.tsv \
    --overlap 0.8 \
    --pred-format minimap
```

### MUMmer4
#### run MUMmer4 on the synthesized data:
```bash
nucmer -p out data/synthetic/synthetic_genome.fasta data/synthetic/synthetic_genome.fasta
mummerplot --fat --png -p out out.delta
show-coords -rcl out.delta > out.coords
show-coords -rcl -T out.delta > results/synthetic/mummer_syn_out.tsv
```
#### MUMmer4 accuracy
```bash
python scripts/evaluate_repeats.py \
    --truth data/synthetic/ground_truth_repeats.tsv \
    --pred results/synthetic/mummer_syn_out.tsv \
    --truth-format groundtruth \
    --pred-format mummer \
    --overlap 0.8
```

## 4. find the long repeats in the large data: the Arabidopsis_thaliana chromosome

To save time and resources, only one chromosome is used to find long repeats. Chromosome 1 of Arabidopsis_thaliana is used in this case.

```bash
python scripts/repeatfinder.py find \
    --genome data/arabidopsis/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa \
    --out results/arabidopsis/arabidopsis_chr1_step4.tsv \
    --chromosome 1 \
    --mask 11110111101111 \
    --step 4 \
    --nseeds 15 \
    --max-occ 300 \
    --band 400 \
    --minlen 20000
```
Similarly, this would also generate a `.tsv` containing long repeats found.

## 5. visualization of Arabidopsis_thaliana results

Since Arabidopsis_thaliana does not have a ground truth of long repeat in the genome, we provide visualization methods based on the results:
```bash
python scripts/visualize_for_our_software.py \
    results/arabidopsis/arabidopsis_chr1_step4.tsv \
    --out results/arabidopsis/plots/our_alg
```

## 6. Benchmark of minimap & MUMmer4 in Arabidopsis_thaliana chromosome

To further inspect how **Minimap** and **MUMmer4** perform in Arabidopsis_thaliana chromosome, we also did the analysis and visualization respectively.

### Minimap2

```bash
samtools faidx data/arabidopsis/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa 1 > data/arabidopsis/chr1.fa

minimap2 -x asm5 -t 1 \
    -N 1000 \
    -p 0.0 \
    data/arabidopsis/chr1.fa data/arabidopsis/chr1.fa \
    > results/arabidopsis/chr1.self.minimap.paf
```
After **Minimap2** out put the result in `.paf` file, we can visualize the result:
```bash
python scripts/visualize_for_minimap.py \
    results/arabidopsis/chr1.self.minimap.paf \
    --out results/arabidopsis/plots/minimap
```
This command would generate a `_transform.tsv` in the current directory, convenient for further analysis or other visualization. Also, 4 plots of the results is generated.

### MUMmer4

```bash
samtools faidx data/arabidopsis/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa 1 > data/arabidopsis/chr1.fa

nucmer -p out data/arabidopsis/chr1.fa data/arabidopsis/chr1.fa
mummerplot --fat --png -p out out.delta
show-coords -rcl out.delta > out.coords
show-coords -rcl -T out.delta > results/arabidopsis/mummer_chr1_out.tsv
```
After **MUMmer4** generated `.csv` as the results, we can visualize them:
```bash
python scripts/visualize_for_mummer.py \
    results/arabidopsis/mummer_chr1_out.tsv \
    --out results/arabidopsis/plots/mummer
```
This command generates 4 plots.
