import random
import argparse
from pathlib import Path

# Fixed seed ensures reproducible synthetic genomes
random.seed(42)


# ----------------------------------------------------------------------
# FASTA utilities
# ----------------------------------------------------------------------
def wrap_fasta(seq: str, width: int = 60) -> str:
    """Wrap DNA sequence into fixed-width FASTA lines."""
    return "\n".join(seq[i:i+width] for i in range(0, len(seq), width))


# ----------------------------------------------------------------------
# Genome generation and mutation models
# ----------------------------------------------------------------------
def generate_random_genome(length: int, gc: float = 0.42) -> str:
    """Generate a random genome with specified GC content."""
    gcb = ['G', 'C']
    atb = ['A', 'T']
    return "".join(
        random.choice(gcb if random.random() < gc else atb)
        for _ in range(length)
    )


def mutate_with_subs_and_indels(seq: str,
                                sub_rate: float = 0.03,
                                indel_rate: float = 0.002,
                                max_indel: int = 3) -> str:
    """Apply substitutions and small indels to a DNA sequence."""
    bases = ['A', 'C', 'G', 'T']
    out = []
    i = 0

    while i < len(seq):

        # indel event
        if random.random() < indel_rate:
            k = random.randint(1, max_indel)
            if random.random() < 0.5 and i + k < len(seq):
                i += k
                continue
            else:
                out.append("".join(random.choice(bases) for _ in range(k)))

        # substitution / no-change
        b = seq[i]
        if random.random() < sub_rate:
            out.append(random.choice([x for x in bases if x != b]))
        else:
            out.append(b)

        i += 1

    return "".join(out)


def revcomp(seq: str) -> str:
    """Reverse complement of a DNA sequence."""
    comp = str.maketrans("ACGT", "TGCA")
    return seq.translate(comp)[::-1]


def insert_fragment(genome: str, pos: int, frag: str):
    """Insert fragment at position pos and return (new_genome, interval)."""
    new_genome = genome[:pos] + frag + genome[pos:]
    return new_genome, (pos, pos + len(frag))


# ----------------------------------------------------------------------
# Dataset builder
# ----------------------------------------------------------------------
def build_dataset(base_len=2_000_000,
                  num_repeats=20,
                  sub_rate=0.03,
                  indel_rate=0.002,
                  allow_reverse=True):

    genome = generate_random_genome(base_len)
    used = []
    records = []

    for rid in range(1, num_repeats + 1):

        # Pick repeat length
        L = random.randint(15_000, 100_000)

        # Find non-overlapping source region
        while True:
            src = random.randint(0, base_len - L - 100_000)
            overlap = any(not (src + L < u0 or src > u1) for u0, u1 in used)
            if not overlap:
                break

        tgt = random.randint(src + L + 50_000, base_len - L - 10_000)

        source = genome[src:src+L]

        # Small random perturbation on mutation rates
        s_rate = max(0, sub_rate + random.uniform(-0.01, 0.01))
        i_rate = max(0, indel_rate + random.uniform(-0.001, 0.001))
        mutated = mutate_with_subs_and_indels(source, s_rate, i_rate,
                                              max_indel=random.randint(2, 5))

        strand = "+"
        if allow_reverse and random.random() < 0.3:
            mutated = revcomp(mutated)
            strand = "-"

        genome, (t0, t1) = insert_fragment(genome, tgt, mutated)

        used.append((src, src + L))
        used.append((t0, t1))

        records.append(dict(
            repeat_id=rid,
            source_start_1based=src + 1,
            source_end_1based=src + L,
            target_start_1based=t0 + 1,
            target_end_1based=t1,
            strand=strand,
            source_length=L,
            inserted_length=len(mutated),
            sub_rate=round(s_rate, 4),
            indel_rate=round(i_rate, 5),
        ))

    return genome, records


# ----------------------------------------------------------------------
# Output writing
# ----------------------------------------------------------------------
def write_outputs(outdir: Path, genome: str, records: list):
    outdir.mkdir(parents=True, exist_ok=True)

    # 1. Complete genome
    with open(outdir / "synthetic_genome.fasta", "w") as f:
        f.write(">chr1\n")
        f.write(wrap_fasta(genome))

    # 2. Ground truth repeats
    cols = ["repeat_id", "source_start_1based", "source_end_1based",
            "target_start_1based", "target_end_1based", "strand",
            "source_length", "inserted_length", "sub_rate", "indel_rate"]

    with open(outdir / "ground_truth_repeats.tsv", "w") as f:
        f.write("\t".join(cols) + "\n")
        for r in records:
            f.write("\t".join(str(r[c]) for c in cols) + "\n")

    # 3. Fragments
    fragdir = outdir / "fragments"
    fragdir.mkdir(exist_ok=True)
    gl = len(genome)
    frag_size = 100_000
    positions = [
        0,
        max(0, gl // 4 - frag_size // 2),
        max(0, gl // 2 - frag_size // 2),
        max(0, 3 * gl // 4 - frag_size // 2),
        max(0, gl - frag_size),
    ]

    for i, s in enumerate(positions, 1):
        e = min(gl, s + frag_size)
        with open(fragdir / f"fragment_{i}.fasta", "w") as f:
            f.write(f">chr1_fragment_{i}:{s+1}-{e}\n")
            f.write(wrap_fasta(genome[s:e]))


# ----------------------------------------------------------------------
# CLI
# ----------------------------------------------------------------------
def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate synthetic genome + ground truth repeat table."
    )
    parser.add_argument("--out", type=str, required=True,
                        help="Output directory (e.g. /data/synthetic/)")
    return parser.parse_args()


# ----------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------
if __name__ == "__main__":
    args = parse_args()
    outdir = Path(args.out)

    genome, records = build_dataset()
    write_outputs(outdir, genome, records)

    print(f"[INFO] Synthetic genome written to:       {outdir / 'synthetic_genome.fasta'}")
    print(f"[INFO] Ground truth written to:          {outdir / 'ground_truth_repeats.tsv'}")
    print(f"[INFO] Fragments stored under:           {outdir / 'fragments'}")