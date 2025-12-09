# BMI3_miniproject

FindLongRepeat is a code developmed by our group to find long repeat in a certain genome.

The whole pipeline is as follows:

## 1. generate synthesized data

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

## 3. accuracy

`evaluate_repeats.py` is used to evaluate the results of the generated long repeats, comparing to the ground truth synthesized.
```bash
python evaluate_repeats.py \
    --truth ground_truth_repeats.tsv \
    --pred results/far-fast.tsv \
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

```
