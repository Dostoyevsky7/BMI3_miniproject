# BMI3_miniproject

FindLongRepeat is a code developmed by our group to find long repeat in a certain genome.

The whole pipeline is as follows:

## 1. generate synthesized data

First, we used data.py to generate a small test dataset by running:
```bash
python data.py
```
This script creates a synthetic genome containing 20 implanted repeats (synthetic_genome.fasta) along with the corresponding ground-truth list (ground_truth_repeats.tsv).

## 2. find long repeat in the synthesized data

Based on the synthesized data, the core code `repeatfinder.py` could be executed to find long repeat:
```bash
python repeatfinder.py find \
    --genome synthetic_genome.fasta \
    --out results/far-fast.tsv \
    --chromosome 1 \
    --mask 11110111101111 \
    --step 4 \
    --nseeds 15 \
    --max-occ 300 \
    --band 400 \
    --minlen 20000
```
This execution would generate a .tsv file that contains the long repeat results.
This script identifies repeat regions within the input genome. Given a dataset as input, it produces a TSV file containing the coordinates, lengths and number of repeats, together with a homology score for each repeat, as well as a summary report (TXT) describing overall results.
When running the tool, the user must specify the input file and a set of parameters. These parameters govern the speed, sensitivity, and memory usage of the algorithm. To make the tool more convenient, we provide several predefined parameter presets—such as fast and balanced. Users may either rely on these presets or manually set their own parameter combinations.

## 3. accuracy

`evaluate_repeats.py` is used to evaluate the results of the generated long repeats, comparing to the ground truth synthesized.
```bash
python evaluate_repeats.py \
    --truth ground_truth_repeats.tsv \
    --pred results/far-fast.tsv \
    --overlap 0.8
```

## 4. Benchmark: other package performance

**minimap2** and **mummer4** are 2 tools rely on **seed-based approach**:

- Minimap2 is a fast, flexible long-read and genome aligner (https://github.com/lh3/minimap2)
- MUMmer4 is a high-precision genome comparison toolkit (https://github.com/mummer4/mummer)

When aligning a genome to itself, any repeated segment will appear as an off-diagonal match, making both tools powerful for detecting long genomic repeats.

Since our algorithm is also **seed-based**, we compare our algorthm with minimap2 and mummer4.

### minimap2
#### run minimap2 on the synthesized data:
```bash
minimap2 -x asm5 -t 1 \
-N 1000 \
-p 0.0   \
synthetic_genome.fasta synthetic_genome.fasta > syn.self.asm5.paf
```
#### minimap2 accuracy
**minimap2** generates a `.paf` file as the result. For the convenience to calculate the accuracy, following command could transform the result into a `.tsv` file:
```bash
python transform.py -i syn.self.asm5.paf -o syn_transform.tsv
```
After transformed into `.tsv`, the accuracy could be calculated by:
```bash
python evaluate_repeats.py \
    --truth ground_truth_repeats.tsv\
    --pred syn_transform.tsv\
    --overlap 0.8\
    --pred-format minimap
```

### mummer4
#### run mummer4 on the synthesized data:
```bash
nucmer -p out synthetic_genome.fasta synthetic_genome.fasta
mummerplot --fat --png -p out out.delta
show-coords -rcl out.delta > out.coords
show-coords -rcl -T out.delta > syn_out.tsv
```
#### mummer4 accuracy
```bash
python evaluate_repeats.py \
    --truth ground_truth_repeats.tsv \
    --pred syn_out.tsv \
    --truth-format groundtruth \
    --pred-format mummer \
    --overlap 0.8
```

## 4. find the long repeats in the large data: the Arabidopsis_thaliana chromosome

To save time and resources, only one chromosome is used to find long repeats. Chromosome 1 of Arabidopsis_thaliana is used in this case.

```bash
python repeatfinder.py find \
    --genome Arabidopsis_thaliana.TAIR10.dna.toplevel.fa \
    --out results/arabidopsis_chr1_step4.tsv \
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
python visualize_for_our_software.py results/arabidopsis_chr1_step4.tsv
```

## 6. Benchmark of minimap & mummer in Arabidopsis_thaliana chromosome

To further inspect how **minimap** and **mummer** perform in Arabidopsis_thaliana chromosome, we also did the analysis and visualization respectively.

### minimap2

```bash
samtools faidx Arabidopsis_thaliana.TAIR10.dna.toplevel.fa 1 > chr1.fa
minimap2 -x asm5 -t 1 \
-N 1000 \
-p 0.0   \
chr1.fa chr1.fa > chr1.self.asm5.paf
```
After **minimap2** out put the result in `.paf` file, we can visualize the result:
```bash
python visualize_for_minimap.py chr1.self.asm5.paf
```
This command would generate a `_transform.tsv` in the current directory, convenient for further analysis or other visualization. Also, 4 plots of the results is generated.

### mummer4

```bash
samtools faidx Arabidopsis_thaliana.TAIR10.dna.toplevel.fa 1 > chr1.fa
nucmer -p out chr1.fa chr1.fa
mummerplot --fat --png -p out out.delta
show-coords -rcl out.delta > out.coords
show-coords -rcl -T out.delta > chr1_out.tsv
```
After **mummer4** generated `.csv` as the results, we can visualize them:
```bash
python visualize_for_mummer.py chr1_out.tsv
```
This command generates 4 plots.
