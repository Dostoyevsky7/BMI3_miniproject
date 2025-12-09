import random
from pathlib import Path

# Set random seed for reproducibility across different runs
# This ensures that the same "random" genome is generated every time
random.seed(42)

def wrap_fasta(seq: str, width: int = 60) -> str:
    """
    Wrap a long DNA sequence into FASTA format with fixed line width.
    
    FASTA files conventionally have sequences wrapped at 60-80 characters per line
    for readability and compatibility with various bioinformatics tools.
    
    Args:
        seq: DNA sequence string (can be any length)
        width: Number of characters per line (default: 60, standard FASTA format)
        
    Returns:
        String with newlines inserted every 'width' characters
        
    Example:
        wrap_fasta("ACGTACGTACGT", width=4) returns "ACGT\\nACGT\\nACGT"
    """
    return "\n".join(seq[i:i+width] for i in range(0, len(seq), width))

def generate_random_genome(length: int, gc: float = 0.42) -> str:
    """
    Generate a random DNA sequence with specified GC content.
    
    GC content (proportion of G and C bases) affects genome properties like
    melting temperature, gene density, and repeat detectability. Real genomes
    vary: humans ~41%, Arabidopsis ~36%, some bacteria >70%.
    
    Args:
        length: Total number of base pairs to generate
        gc: Target GC content as a fraction (0.0 to 1.0)
            Default 0.42 approximates Arabidopsis thaliana
            
    Returns:
        Random DNA sequence string containing only A, C, G, T
        
    Algorithm:
        For each position, randomly decide if it should be GC or AT
        based on gc parameter, then randomly pick one base from that pair
        
    Example:
        generate_random_genome(100, gc=0.5) generates 100bp with ~50% GC
    """
    gcb = ['G', 'C']  # GC base pair
    atb = ['A', 'T']  # AT base pair
    return "".join(
        random.choice(gcb if random.random() < gc else atb) 
        for _ in range(length)
    )

def mutate_with_subs_and_indels(
    seq: str, 
    sub_rate: float = 0.03, 
    indel_rate: float = 0.002, 
    max_indel: int = 3
) -> str:
    """
    Simulate evolutionary mutations on a DNA sequence.
    
    This function models realistic genomic mutations that occur during DNA replication
    and over evolutionary time:
    - Substitutions: Single base changes (e.g., A->G)
    - Insertions: Addition of 1-N bases
    - Deletions: Removal of 1-N bases
    
    Args:
        seq: Input DNA sequence to mutate
        sub_rate: Probability of substitution per base (default: 3% = 0.03)
                 Realistic values: 0.01-0.05 for diverged duplicates
        indel_rate: Probability of insertion/deletion per base (default: 0.2% = 0.002)
                   Indels are less common than substitutions in real genomes
        max_indel: Maximum length of a single indel event (default: 3bp)
                  Small indels (1-5bp) are most common in nature
                  
    Returns:
        Mutated DNA sequence (length may differ from input due to indels)
        
    Algorithm:
        1. Iterate through each base in input sequence
        2. At each position, randomly decide:
           a. Apply indel (with probability indel_rate)
              - 50% chance: deletion (skip 1-max_indel bases)
              - 50% chance: insertion (add 1-max_indel random bases)
           b. Apply substitution (with probability sub_rate)
              - Replace current base with different base
           c. Keep original base (1 - sub_rate - indel_rate)
        3. Return concatenated result
        
    Biological rationale:
        - Substitutions: point mutations from replication errors, UV damage
        - Insertions/deletions: slippage during replication, transposon activity
        - Indel rate << sub rate: DNA polymerase is more accurate at maintaining length
        
    Example:
        Input:  "ACGTACGT" (8bp)
        Output: "ACGAACGT" (8bp, one substitution T->A)
        Output: "ACGACGT"  (7bp, one deletion of T)
        Output: "ACGTTACGT" (9bp, one insertion of T)
    """
    bases = ['A', 'C', 'G', 'T']
    out = []
    i = 0
    while i < len(seq):
        # Randomly apply indel (insertion or deletion)
        if random.random() < indel_rate:
            k = random.randint(1, max_indel)  # Length of indel
            if random.random() < 0.5 and i + k < len(seq):  
                # Deletion: skip next k bases
                i += k
                continue
            else:  
                # Insertion: add k random bases
                out.append("".join(random.choice(bases) for _ in range(k)))
        
        # Process current base
        b = seq[i]
        if random.random() < sub_rate:
            # Substitution: replace with different base
            out.append(random.choice([x for x in bases if x != b]))
        else:
            # No mutation: keep original base
            out.append(b)
        i += 1
    
    return "".join(out)

def revcomp(seq: str) -> str:
    """
    Compute the reverse complement of a DNA sequence.
    
    In double-stranded DNA, the two strands are antiparallel and complementary:
    - A pairs with T
    - C pairs with G
    
    Reverse complement is essential for:
    - Detecting repeats on opposite strands
    - Primer design
    - Understanding gene orientation
    
    Args:
        seq: DNA sequence string (5' to 3' direction)
        
    Returns:
        Reverse complement sequence (3' to 5' direction, then reversed to 5' to 3')
        
    Algorithm:
        1. Complement: A↔T, C↔G (using str.translate)
        2. Reverse: flip the string (using [::-1])
        
    Example:
        revcomp("ACGT") returns "ACGT" (palindrome)
        revcomp("AAAA") returns "TTTT"
        revcomp("ACGT") = "ACGT" (complement: TGCA, reverse: ACGT)
        
    Note: For real genomes, also handle ambiguous bases (N, R, Y, etc.)
    """
    comp = str.maketrans("ACGT", "TGCA")  # Complementarity table
    return seq.translate(comp)[::-1]      # Complement then reverse

def insert_fragment(genome: str, start0: int, fragment: str) -> (str, tuple):
    """
    Insert a DNA fragment into a genome at a specified position (WITHOUT replacement).
    
    This simulates tandem duplication or transposon insertion events where
    new sequence is added rather than replacing existing sequence. This causes
    the genome to grow in length.
    
    Args:
        genome: Original genome sequence
        start0: Insertion position (0-based index)
               Fragment will be inserted BEFORE this position
        fragment: DNA sequence to insert
        
    Returns:
        Tuple of (new_genome, inserted_interval):
        - new_genome: Modified genome with fragment inserted
        - inserted_interval: (start, end) in 0-based half-open coordinates
                            where fragment now resides in new_genome
                            
    Coordinate system:
        - Input: 0-based (Python standard)
        - Output interval: [start, end) half-open (standard in genomics)
        
    Example:
        genome = "AAAA", fragment = "GG", start0 = 2
        Result: ("AAGGAA", (2, 4))
        Explanation: "AA" + "GG" + "AA", fragment occupies positions [2,4)
        
    Biological context:
        Real duplications can occur via:
        - Tandem duplication (insertion near source)
        - Transposon-mediated duplication (insertion anywhere)
        - Unequal crossing over during meiosis
    """
    new_genome = genome[:start0] + fragment + genome[start0:]
    # Fragment now occupies [start0, start0+len(fragment))
    return new_genome, (start0, start0 + len(fragment))

def build_dataset(
    base_len=2_000_000,   # Base genome length before insertions (2 Mbp)
    dup_lengths=None,      # Legacy parameter (now unused, kept for compatibility)
    sub_rate=0.03,        # Default substitution rate (3% divergence)
    indel_rate=0.002,     # Default indel rate (0.2% frequency)
    allow_reverse=True,   # Allow reverse-complemented duplications
    num_repeats=20        # Number of repeat regions to insert
):
    """
    Generate a synthetic genome with realistic tandem duplications for testing.
    
    This function creates a controlled test dataset where:
    - The "true" repeat locations are known (ground truth)
    - Repeat properties (length, divergence, orientation) are diverse
    - The genome is large enough to be realistic but small enough to process quickly
    
    Typical use case:
        Test a repeat-finding algorithm by comparing its output to ground_truth_repeats.tsv
        
    Args:
        base_len: Initial genome size in base pairs (default: 2 Mbp)
                 Realistic for bacterial genomes or eukaryotic chromosomes
                 Final size will be larger due to insertions
                 
        dup_lengths: Deprecated parameter (kept for backward compatibility)
                    Repeat lengths are now generated dynamically
                    
        sub_rate: Base substitution rate for mutations (default: 0.03 = 3%)
                 Will be randomly perturbed ±1% per repeat to add diversity
                 
        indel_rate: Base insertion/deletion rate (default: 0.002 = 0.2%)
                   Will be randomly perturbed ±0.1% per repeat
                   
        allow_reverse: If True, 30% of repeats will be reverse-complemented
                      Simulates duplications followed by inversion events
                      
        num_repeats: Number of duplicate regions to insert (default: 20)
                    More repeats = more realistic but slower processing
                    
    Returns:
        Tuple of (genome, records):
        - genome: Final DNA sequence string (length > base_len)
        - records: List of dictionaries, one per repeat, containing:
            * repeat_id: Unique identifier (1-indexed)
            * source_start_1based: Start of source region (1-based, inclusive)
            * source_end_1based: End of source region (1-based, inclusive)
            * target_start_1based: Start of inserted copy (1-based, inclusive)
            * target_end_1based: End of inserted copy (1-based, inclusive)
            * strand: "+" (forward) or "-" (reverse complement)
            * source_length: Length of source region in bp
            * inserted_length: Length of inserted copy (may differ due to indels)
            * sub_rate: Actual substitution rate used for this repeat
            * indel_rate: Actual indel rate used for this repeat
            
    Algorithm:
        1. Generate random base genome (no repeats)
        2. For each of num_repeats iterations:
           a. Randomly choose repeat length (15-100 kb)
           b. Randomly choose source position (avoid overlaps with existing repeats)
           c. Randomly choose target position (after source, with buffer)
           d. Extract source sequence
           e. Mutate copy with randomized rates (adds diversity)
           f. Optionally reverse-complement (30% chance)
           g. Insert copy into genome at target
           h. Record metadata for ground truth
        3. Return final genome and all records
        
    Design considerations:
        - Overlap prevention: Repeats don't overlap, making ground truth unambiguous
        - Length diversity: 15-100 kb covers typical gene-scale to operon-scale
        - Divergence diversity: Random perturbations prevent all repeats from
          having identical identity scores
        - Spatial distribution: Target always after source with buffer prevents
          trivial self-matches and creates realistic tandem duplications
        - Strand diversity: 30% reverse matches real inversion frequency in genomes
        
    Validation uses:
        - Precision: How many predicted repeats are in ground truth?
        - Recall: How many ground truth repeats were found?
        - Boundary accuracy: How close are predicted coordinates to truth?
        
    Example output:
        genome length: 2,400,000 bp (grew by 400kb due to 20 insertions)
        records: 20 entries with diverse properties
    """
    # Step 1: Generate base genome with no repeats
    genome = generate_random_genome(base_len)
    records = []
    
    # Track used regions to prevent overlaps (ensures unambiguous ground truth)
    used_regions = []  # List of (start, end) tuples
    
    # Step 2: Insert num_repeats duplicate regions
    for i in range(1, num_repeats + 1):
        # Randomly generate repeat length (15kb - 100kb range)
        # This covers typical gene clusters, operons, and large tandem duplications
        L = random.randint(15_000, 100_000)
        
        # Find non-overlapping source position
        # Retry until we find a region that doesn't overlap with existing repeats
        while True:
            # Leave 100kb buffer at end for target insertion
            src = random.randint(0, base_len - L - 100_000)
            
            # Check for overlap with any existing used region
            overlap = any(
                not (src + L < r[0] or src > r[1])  # Intervals overlap if neither is completely before the other
                for r in used_regions
            )
            if not overlap:
                break  # Found a valid non-overlapping source
        
        # Choose target position for the duplicate
        # Must be:
        # - After source (src + L + 50kb buffer)
        # - Before genome end (leave 10kb buffer)
        # This creates realistic "tandem duplication" geometry
        tgt = random.randint(src + L + 50_000, base_len - L - 10_000)
        
        # Extract source sequence to be duplicated
        source_seq = genome[src:src+L]
        
        # Add evolutionary mutations to the copy (simulate divergence over time)
        # Randomize rates ±1% (sub) and ±0.1% (indel) to add diversity
        # This prevents all repeats from having identical divergence
        actual_sub = sub_rate + random.uniform(-0.01, 0.01)
        actual_indel = indel_rate + random.uniform(-0.001, 0.001)
        copy = mutate_with_subs_and_indels(
            source_seq, 
            sub_rate=max(0, actual_sub),        # Ensure non-negative
            indel_rate=max(0, actual_indel),    # Ensure non-negative
            max_indel=random.randint(2, 5)      # Vary indel size (2-5bp)
        )
        
        # Randomly decide orientation (30% probability of inversion)
        # This simulates duplications followed by inversion events
        strand = "+"
        if allow_reverse and random.random() < 0.3:
            copy = revcomp(copy)
            strand = "-"
        
        # Insert the mutated copy into the genome
        # This grows the genome length by len(copy)
        genome, ins_iv = insert_fragment(genome, tgt, copy)
        
        # Update used regions to include both source and target
        # This prevents future repeats from overlapping with these
        used_regions.append((src, src + L))           # Source region
        used_regions.append((ins_iv[0], ins_iv[1]))   # Inserted copy region
        
        # Record ground truth metadata for evaluation
        # All coordinates are 1-based (biological convention) for output
        rec = dict(
            repeat_id=i,
            source_start_1based=src+1,              # Convert 0-based to 1-based
            source_end_1based=src+L,                # Already inclusive
            target_start_1based=ins_iv[0]+1,        # Convert 0-based to 1-based
            target_end_1based=ins_iv[1],            # Already inclusive
            strand=strand,
            source_length=L,
            inserted_length=len(copy),              # May differ from L due to indels
            sub_rate=round(actual_sub, 4),          # Record actual rate used
            indel_rate=round(actual_indel, 5)       # Record actual rate used
        )
        records.append(rec)
    
    return genome, records

if __name__ == "__main__":
    # Generate synthetic genome and ground truth
    genome, records = build_dataset()

    # Set output directory (current directory)
    out_dir = Path(".")
    
    # Create subdirectory for fragment files
    (out_dir / "fragments").mkdir(parents=True, exist_ok=True)

    # 1) Write complete genome to FASTA file
    #    Standard FASTA format: >header followed by wrapped sequence
    with open(out_dir / "synthetic_genome.fasta", "w") as f:
        f.write(">chr1\n")  # Single chromosome named "chr1"
        f.write(wrap_fasta(genome, 60))  # Wrap at 60 characters per line

    # 2) Write ground truth repeat coordinates to TSV file
    #    This file is used for evaluating repeat-finding algorithm accuracy
    with open(out_dir / "ground_truth_repeats.tsv", "w") as f:
        # Define column order
        cols = ["repeat_id","source_start_1based","source_end_1based",
                "target_start_1based","target_end_1based","strand",
                "source_length","inserted_length","sub_rate","indel_rate"]
        
        # Write header
        f.write("\t".join(cols) + "\n")
        
        # Write data rows (one per repeat)
        for r in records:
            f.write("\t".join(str(r[c]) for c in cols) + "\n")

    # 3) Write 5 genome fragments for quick parameter testing
    #    Rationale: Testing on small fragments (100kb each) is much faster
    #    than processing the entire 2+ Mbp genome. Useful for:
    #    - Quick parameter optimization (e.g., testing different seed lengths)
    #    - Debugging algorithm logic
    #    - Iterative development
    #
    #    Fragment selection strategy:
    #    - Fragment 1: Start of genome (position 0)
    #    - Fragment 2: 1/4 through genome
    #    - Fragment 3: Middle of genome
    #    - Fragment 4: 3/4 through genome
    #    - Fragment 5: End of genome
    #
    #    This ensures fragments sample different parts of the genome,
    #    likely capturing different repeat types and densities.
    
    gl = len(genome)
    frag_size = 100_000  # 100kb per fragment (manageable size)
    
    # Calculate start positions for 5 evenly-spaced fragments
    starts = [
        0,                                    # Beginning
        max(0, gl//4 - frag_size//2),        # 1/4 point (centered)
        max(0, gl//2 - frag_size//2),        # 1/2 point (centered)
        max(0, 3*gl//4 - frag_size//2),      # 3/4 point (centered)
        max(0, gl - frag_size)               # End (last 100kb)
    ]
    
    # Write each fragment to a separate FASTA file
    for idx, s in enumerate(starts, 1):
        e = min(s + frag_size, gl)  # Ensure we don't exceed genome length
        
        with open(out_dir / "fragments" / f"fragment_{idx}.fasta", "w") as f:
            # Header includes genomic coordinates for traceability
            # Format: >chr1_fragment_N:start-end (1-based, inclusive)
            f.write(f">chr1_fragment_{idx}:{s+1}-{e}\n")
            f.write(wrap_fasta(genome[s:e], 60))
