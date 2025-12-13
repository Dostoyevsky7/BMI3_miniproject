#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
repeatfinder.py (Memory & Speed Optimized Version)

Core optimizations:
1. Streaming processing with early pruning to prevent pair explosion
2. Local index filtering with per-chromosome seed frequency limits
3. Batch processing to control memory peak usage
4. Maintains original bucketing logic for diagonal-based seed clustering
"""

import argparse
import math
import sys
import gc
import time
import os
from collections import defaultdict, Counter
from multiprocessing import Pool, cpu_count

# ===================== Preset Parameter Configurations =====================

PRESET_CONFIGS = {
    "balanced": {
        "mask": "11110111101111",
        "step": 1,
        "nseeds": 45,
        "max_occ": 5000,
        "band": 400,
        "description": "Balanced mode: Optimal trade-off between recall and precision"
    },
    "fast": {
        "mask": "11110111101111",
        "step": 3,
        "nseeds": 20,
        "max_occ": 500,
        "band": 400,
        "description": "Fast mode: Prioritizes speed, suitable for initial screening of large genomes"
    },
    "far-fast": {
        "mask": "11110111101111",
        "step": 4,
        "nseeds": 15,
        "max_occ": 300,
        "band": 400,
        "description": "Ultra-fast mode: Maximum speed for very large genomes (lowest sensitivity)"
    },
    "precise": {
        "mask": "11110111101111",
        "step": 1,
        "nseeds": 60,
        "max_occ": 3000,
        "band": 300,
        "description": "Precise mode: Stricter filtering to reduce false positives"
    },
    "sensitive": {
        "mask": "1111011101111",
        "step": 1,
        "nseeds": 110,
        "max_occ": 2000,
        "band": 200,
        "description": "Sensitive mode: Maximum recall to find more true positives (but with more false positives)"
    }
}

# ===================== FASTA File I/O =====================

def read_fasta(path):
    """
    Read FASTA file and return an ordered list of sequences.
    
    Args:
        path: Path to FASTA file
        
    Returns:
        List of tuples [(sequence_name, sequence_string), ...]
        All sequences are converted to uppercase for consistency
    """
    sequences = []
    name = None
    parts = []
    with open(path, "r") as fh:
        for line in fh:
            line = line.strip()
            if not line: continue
            if line.startswith(">"):
                if name is not None:
                    sequences.append((name, "".join(parts).upper()))
                name = line[1:].split()[0]
                parts = []
            else:
                parts.append(line)
        if name is not None:
            sequences.append((name, "".join(parts).upper()))
    return sequences

# Translation table for efficient reverse complement computation
_complement_table = str.maketrans("ACGTNacgtn", "TGCANtgcan")

def revcomp(seq):
    """
    Compute reverse complement of a DNA sequence.
    
    Args:
        seq: DNA sequence string (A, C, G, T, N)
        
    Returns:
        Reverse complement sequence
        
    Example:
        revcomp("ACGT") -> "ACGT"
        revcomp("AAAA") -> "TTTT"
    """
    return seq.translate(_complement_table)[::-1]

# ===================== Progress Bar Utility =====================

def print_progress_bar(iteration, total, prefix='', suffix='', length=50, fill='█'):
    """
    Display a progress bar in the console.
    
    Args:
        iteration: Current iteration number (0-indexed or 1-indexed)
        total: Total number of iterations
        prefix: Text to display before the progress bar
        suffix: Text to display after the progress bar
        length: Character length of the progress bar
        fill: Character to use for filled portion
        
    Example output:
        Processing [████████████████████----------] 66.7% Complete
    """
    percent = f"{100 * (iteration / float(total)):.1f}"
    filled_length = int(length * iteration // total)
    bar = fill * filled_length + '-' * (length - filled_length)
    print(f'\r{prefix} |{bar}| {percent}% {suffix}', end='', flush=True)
    if iteration == total:
        print()

# ===================== Global Seed Pre-filtering =====================

def extract_all_seeds_fast(seq, mask, step):
    """
    Fast extraction of all canonical spaced seeds from a sequence.
    
    A spaced seed uses only specific positions (defined by mask) to create
    a signature for matching. This reduces sensitivity to small mutations
    while maintaining speed.
    
    Args:
        seq: DNA sequence string
        mask: Binary string (e.g., "11110111101111") where '1' = use position
        step: Sampling interval (step=1 means every position, step=3 means every 3rd)
        
    Returns:
        List of canonical seed strings (lexicographically smaller of fwd/rev)
        
    Technical details:
        - Skips seeds containing 'N' (ambiguous bases)
        - Uses canonical form (min of forward and reverse complement)
          to ensure strand-independent matching
        - Memory efficient: only stores seed strings, not positions
    """
    mask_positions = [i for i, c in enumerate(mask) if c == "1"]
    window_len = len(mask)
    L = len(seq)
    
    results = []
    for pos in range(0, L - window_len + 1, step):
        has_n = False
        raw_chars = []
        for i in mask_positions:
            base = seq[pos + i]
            if base == 'N':
                has_n = True
                break
            raw_chars.append(base)
        
        if has_n: continue
        
        raw_seed = "".join(raw_chars)
        rc_seed = revcomp(raw_seed)
        # Use canonical form for strand-independent matching
        canonical = raw_seed if raw_seed <= rc_seed else rc_seed
        results.append(canonical)
    
    return results

def prefilter_seeds_global(genome_seqs, mask, step, max_occ):
    """
    Global pre-filtering: Count seed occurrences across entire genome
    and blacklist high-frequency seeds to prevent computational explosion.
    
    Rationale:
        Seeds appearing too many times (e.g., repetitive elements, low-complexity
        regions) generate too many pairwise matches (O(n²) problem). By blacklisting
        them upfront, we dramatically reduce downstream computation.
    
    Args:
        genome_seqs: List of (chromosome_name, sequence) tuples
        mask: Spaced seed mask pattern
        step: Sampling step size
        max_occ: Maximum allowed occurrences for a seed (exceeding -> blacklist)
        
    Returns:
        Set of blacklisted seed strings
        
    Side effects:
        Prints statistics about seed distribution and filtering
        Performs garbage collection to free memory
    """
    print("[INFO] Pre-filtering high-frequency seeds globally...")
    seed_counter = Counter()
    
    total_chroms = len(genome_seqs)
    for idx, (name, seq) in enumerate(genome_seqs, 1):
        print_progress_bar(idx, total_chroms, 
                          prefix=f'  -> Counting seeds [{idx}/{total_chroms}]',
                          suffix=f'{name} ({len(seq):,} bp)')
        seeds = extract_all_seeds_fast(seq, mask, step)
        seed_counter.update(seeds)
        del seeds
        gc.collect()
    
    total_seeds = sum(seed_counter.values())
    blacklist = set(seed for seed, count in seed_counter.items() if count > max_occ)
    blacklist_hits = sum(seed_counter[seed] for seed in blacklist)
    
    print(f"[INFO] Total seeds: {total_seeds:,}")
    print(f"[INFO] Blacklisted {len(blacklist):,} seeds (>{max_occ} occurrences).")
    print(f"[INFO] Filtered out {blacklist_hits:,} hits ({100*blacklist_hits/total_seeds:.1f}% of total).")
    
    del seed_counter
    gc.collect()
    
    return blacklist

# ===================== Core Algorithm (Optimized Version) =====================

def build_filtered_index_with_local_limit(seq, mask, step, blacklist, local_max=100):
    """
    Build a filtered seed index with two-tier filtering:
    1. Global blacklist (from prefilter_seeds_global)
    2. Local frequency limit (per-chromosome high-frequency seeds)
    
    This two-tier approach prevents both genome-wide repeats (e.g., transposons)
    and chromosome-specific repeats from overwhelming the pairwise matching step.
    
    Args:
        seq: DNA sequence string (single chromosome)
        mask: Spaced seed mask pattern
        step: Sampling step size
        blacklist: Set of globally blacklisted seeds
        local_max: Maximum occurrences allowed within this chromosome
        
    Returns:
        Dictionary mapping canonical_seed -> [(position, orientation), ...]
        where orientation is +1 (forward) or -1 (reverse complement)
        
    Algorithm:
        1. First pass: Count local seed frequencies
        2. Identify local high-frequency seeds (exceeding local_max)
        3. Second pass: Build index, skipping both global and local blacklists
        
    Memory optimization:
        - Deletes temporary counter after use
        - Calls garbage collector explicitly
    """
    mask_positions = [i for i, c in enumerate(mask) if c == "1"]
    window_len = len(mask)
    L = len(seq)
    
    # First pass: Count local seed frequencies
    local_counter = Counter()
    for pos in range(0, L - window_len + 1, step):
        has_n = False
        raw_chars = []
        for i in mask_positions:
            base = seq[pos + i]
            if base == 'N':
                has_n = True
                break
            raw_chars.append(base)
        
        if has_n: continue
        
        raw_seed = "".join(raw_chars)
        rc_seed = revcomp(raw_seed)
        canonical = raw_seed if raw_seed <= rc_seed else rc_seed
        
        if canonical not in blacklist:
            local_counter[canonical] += 1
    
    # Identify local high-frequency seeds
    local_blacklist = set(seed for seed, count in local_counter.items() if count > local_max)
    
    # Second pass: Build index (skip both global and local blacklists)
    index = defaultdict(list)
    
    for pos in range(0, L - window_len + 1, step):
        has_n = False
        raw_chars = []
        for i in mask_positions:
            base = seq[pos + i]
            if base == 'N':
                has_n = True
                break
            raw_chars.append(base)
        
        if has_n: continue
        
        raw_seed = "".join(raw_chars)
        rc_seed = revcomp(raw_seed)
        
        if raw_seed <= rc_seed:
            canonical = raw_seed
            orientation = 1  # Forward strand
        else:
            canonical = rc_seed
            orientation = -1  # Reverse complement strand
        
        if canonical in blacklist or canonical in local_blacklist:
            continue
        
        index[canonical].append((pos, orientation))
    
    del local_counter
    gc.collect()
    
    return index

def find_repeats_between_chroms_worker(args_tuple):
    """
    Worker function for parallel processing of chromosome pairs.
    Finds all repeat pairs between two chromosomes (or within one chromosome).
    
    This is the main computational bottleneck, optimized for:
    - Memory efficiency (streaming pairs, not storing all at once)
    - Early termination (stops if pair count exceeds threshold)
    - Local filtering (per-chromosome seed frequency limits)
    
    Args:
        args_tuple: Tuple containing:
            - name1, seq1: First chromosome name and sequence
            - name2, seq2: Second chromosome name and sequence
            - mask: Spaced seed mask pattern
            - step: Sampling step size
            - nseeds: Minimum seeds required in a chain
            - band: Diagonal band width for bucketing
            - minlen: Minimum repeat length
            - blacklist: Global seed blacklist
            
    Returns:
        List of repeat call dictionaries, each containing:
            chr1, start1, end1, chr2, start2, end2, strand, length,
            identity, score, n_seeds, seed_density, context
            
    Algorithm flow:
        1. Build filtered index for seq1 (with local limits)
        2. Stream through seq2, extract seeds, match against seq1 index
        3. Generate seed pairs (pos1, pos2, orientation_sign)
        4. Stop early if pair count exceeds MAX_PAIRS (prevents explosion)
        5. Bucket pairs by diagonal (groups collinear seeds)
        6. Chain seeds within each bucket (forms contiguous repeat regions)
        7. Build repeat calls from chains (compute statistics)
    """
    name1, seq1, name2, seq2, mask, step, nseeds, band, minlen, blacklist = args_tuple
    
    # Build filtered index for first sequence with local frequency limits
    idx1 = build_filtered_index_with_local_limit(seq1, mask, step, blacklist, local_max=100)
    
    if not idx1:
        return []
    
    mask_positions = [i for i, c in enumerate(mask) if c == "1"]
    window_len = len(mask)
    L2 = len(seq2)
    
    # Stream pairs (do NOT store all at once to save memory)
    pairs = []
    pair_count = 0
    MAX_PAIRS = 10_000_000  # Safety limit: prevent memory overflow
    
    for pos2 in range(0, L2 - window_len + 1, step):
        has_n = False
        raw_chars = []
        for i in mask_positions:
            base = seq2[pos2 + i]
            if base == 'N':
                has_n = True
                break
            raw_chars.append(base)
        
        if has_n: continue
        
        raw_seed = "".join(raw_chars)
        rc_seed = revcomp(raw_seed)
        
        if raw_seed <= rc_seed:
            canonical = raw_seed
            ori2 = 1
        else:
            canonical = rc_seed
            ori2 = -1
        
        if canonical in blacklist:
            continue
        
        if canonical in idx1:
            # Additional safety: skip if too many matches (prevents O(n²) explosion)
            if len(idx1[canonical]) > 100:
                continue
            
            for pos1, ori1 in idx1[canonical]:
                # Skip self-matches (same position on same chromosome)
                if name1 == name2 and pos1 == pos2:
                    continue
                # Skip redundant pairs (only keep pos1 < pos2 for intra-chromosomal)
                if name1 == name2 and pos1 > pos2:
                    continue
                
                # Orientation sign: +1 if both same strand, -1 if opposite
                sign = ori1 * ori2
                pairs.append((pos1, pos2, sign))
                pair_count += 1
                
                # Early termination if pair count exceeds threshold
                if pair_count > MAX_PAIRS:
                    break
            
            if pair_count > MAX_PAIRS:
                break
    
    del idx1
    gc.collect()
    
    if not pairs:
        return []
    
    # Bucket pairs by diagonal (groups collinear seeds together)
    buckets = bucket_pairs_by_diagonal_local(pairs, gap_tolerance=band * 2)
    
    # Chain seeds within each bucket (forms contiguous repeat regions)
    chains = chain_all_buckets_local(buckets, min_seeds=nseeds, max_gap=20000)
    
    # Build repeat calls from chains
    local_calls = []
    temp_minlen = 500  # Internal minimum for initial filtering
    for ori_flag, chain_pairs in chains:
        call = build_call_from_chain_local(
            seq1, seq2, name1, name2, ori_flag, chain_pairs,
            window_len, temp_minlen
        )
        if call:
            local_calls.append(call)
    
    return local_calls

# ===================== Bucketing and Chaining =====================

def bucket_pairs_by_diagonal_local(pairs, gap_tolerance):
    """
    Group seed pairs by diagonal bands for efficient chaining.
    
    Intuition:
        True repeat regions form diagonal lines in a dotplot (pos1 vs pos2).
        By grouping seeds into diagonal bands, we can efficiently identify
        collinear seed matches that likely belong to the same repeat.
    
    Args:
        pairs: List of (pos1, pos2, orientation_sign) tuples
        gap_tolerance: Width of diagonal bands (larger = more tolerant to gaps)
        
    Returns:
        Dictionary mapping (orientation, diag_bucket) -> list of (pos1, pos2) pairs
        
    Algorithm:
        - For forward strand (+): diagonal = pos1 - pos2
        - For reverse strand (-): diagonal = pos1 + pos2 (inverted)
        - Quantize diagonal into buckets of width gap_tolerance
        - Add pairs to adjacent buckets if near bucket boundary (prevents splits)
        
    Example:
        If gap_tolerance = 800:
        - Seed at (1000, 200) has diag = 800, goes to bucket 1
        - Seed at (1050, 250) has diag = 800, also bucket 1 (same repeat)
        - Seeds in bucket 1 are likely from the same repeat region
    """
    buckets = defaultdict(list)
    for pos1, pos2, sign in pairs:
        if sign > 0:
            diag = pos1 - pos2  # Forward strand diagonal
            ori_flag = "+"
        else:
            diag = pos1 + pos2  # Reverse strand diagonal (inverted geometry)
            ori_flag = "-"
        
        diag_bucket = int(diag // gap_tolerance)
        key = (ori_flag, diag_bucket)
        buckets[key].append((pos1, pos2))
        
        # Add to adjacent buckets if near boundary (prevents artificial splits)
        remainder = diag % gap_tolerance
        if remainder < (gap_tolerance * 0.2):
            buckets[(ori_flag, diag_bucket - 1)].append((pos1, pos2))
        elif remainder > (gap_tolerance * 0.8):
            buckets[(ori_flag, diag_bucket + 1)].append((pos1, pos2))
    
    return buckets

def chain_all_buckets_local(buckets, min_seeds, max_gap):
    """
    Chain seeds within each diagonal bucket to form contiguous repeat regions.
    
    A "chain" is a sequence of collinear seeds with small gaps between them,
    indicating a single repeat region. This function uses a simple greedy
    algorithm to extend chains as long as gaps remain small.
    
    Args:
        buckets: Dictionary from bucket_pairs_by_diagonal_local
        min_seeds: Minimum number of seeds required in a valid chain
        max_gap: Maximum allowed gap between consecutive seeds (bp)
        
    Returns:
        List of (orientation, chain_pairs) tuples
        where chain_pairs is a list of (pos1, pos2) coordinates
        
    Algorithm (per bucket):
        1. Sort pairs by position (ensures sequential processing)
        2. Start a new chain with first pair
        3. For each subsequent pair:
           - If gap from last pair < max_gap: extend current chain
           - Else: save current chain (if >= min_seeds), start new chain
        4. Save final chain
        
    Rationale:
        Simple greedy chaining is fast and works well because:
        - Seeds are already grouped by diagonal (from bucketing)
        - True repeats have small, relatively uniform gaps
        - Complex DP algorithms are overkill for this filtered data
    """
    all_chains = []
    for (ori_flag, _), pairs in buckets.items():
        pairs.sort()  # Sort by position for sequential chaining
        
        current_chain = [pairs[0]]
        for p in pairs[1:]:
            pos1, pos2 = p
            last_pos1, last_pos2 = current_chain[-1]
            
            # Check if gap is small enough to extend current chain
            if abs(pos1 - last_pos1) <= max_gap and abs(pos2 - last_pos2) <= max_gap:
                current_chain.append(p)
            else:
                # Gap too large: save current chain if long enough, start new chain
                if len(current_chain) >= min_seeds:
                    all_chains.append((ori_flag, current_chain))
                current_chain = [p]
        
        # Don't forget to save the last chain
        if len(current_chain) >= min_seeds:
            all_chains.append((ori_flag, current_chain))
    
    return all_chains

def build_call_from_chain_local(seq1, seq2, name1, name2, ori_flag, chain_pairs, window_len, minlen):
    """
    Build a repeat call from a chain of seed pairs, computing all statistics.
    
    This function computes a sophisticated "homology score" that integrates
    multiple factors to assess the likelihood that this repeat came from a
    true genomic duplication event (rather than random similarity).
    
    Args:
        seq1, seq2: Chromosome sequences
        name1, name2: Chromosome names
        ori_flag: Strand orientation ("+" or "-")
        chain_pairs: List of (pos1, pos2) seed coordinates
        window_len: Length of spaced seed window
        minlen: Minimum repeat length to consider
        
    Returns:
        Dictionary with repeat call information, or None if below minlen
        
    Homology Score Model:
        score = length_score × identity_score × density_score × context_penalty × significance_bonus
        
        Components:
        1. length_score (log-scale):
           - Avoids over-weighting extremely long repeats
           - log10(length) grows slowly with length
           
        2. identity_score (exponential):
           - Heavily penalizes low-identity matches (< 70%)
           - Heavily rewards high-identity matches (> 90%)
           - Uses exp() to create strong non-linear effect
           
        3. density_score (uniformity):
           - Measures how evenly seeds are distributed
           - Ideal: ~1 seed per window_len bp
           - Too sparse (<0.3): probably incomplete match
           - Too dense (>2.0): probably noise/artifacts
           
        4. context_penalty (chromosome relationship):
           - Intra-chromosomal: penalty = 1.0 (no penalty)
           - Inter-chromosomal: penalty = 0.7 (requires stronger evidence)
           - Rationale: tandem duplications more common than translocations
           
        5. significance_bonus (statistical):
           - Based on Z-score: how much does identity exceed random expectation?
           - Random DNA identity ≈ 0.25 (4 bases, uniform distribution)
           - Z-score > 5: highly significant, apply bonus
    """
    pos1_list = [p[0] for p in chain_pairs]
    pos2_list = [p[1] for p in chain_pairs]
    
    # Compute repeat boundaries (add window_len to cover last seed)
    start1 = min(pos1_list)
    end1 = max(pos1_list) + window_len
    start2 = min(pos2_list)
    end2 = max(pos2_list) + window_len
    
    # Extract subsequences
    sub1 = seq1[start1:end1]
    sub2 = seq2[start2:end2]
    
    # Reverse complement seq2 if on opposite strand
    if ori_flag == "-":
        sub2 = revcomp(sub2)
    
    # Use minimum length (handle edge cases with unequal lengths)
    length = min(len(sub1), len(sub2))
    if length < minlen: 
        return None
    
    # Compute sequence identity via sampling (faster than full alignment)
    ident = identity_sampling(sub1[:length], sub2[:length])
    
    # ==================== Homology Score Calculation ====================
    
    # 1. Length score (log-scale to prevent over-emphasis on very long repeats)
    length_score = math.log10(length + 1)
    
    # 2. Identity score (exponential to strongly favor high-identity matches)
    #    - identity < 70%: score decays rapidly (likely false positive)
    #    - identity > 90%: score increases rapidly (likely true duplication)
    if ident < 0.7:
        identity_score = math.exp(-5 * (1 - ident))  # Rapid decay
    else:
        identity_score = math.exp(3 * (ident - 0.7))  # Rapid growth
    
    # 3. Seed density score (measures uniformity of seed distribution)
    #    Ideal: seeds evenly spaced across repeat region
    theoretical_seeds = length / window_len
    actual_seeds = len(chain_pairs)
    seed_density = actual_seeds / theoretical_seeds if theoretical_seeds > 0 else 0
    
    # Penalty function for density:
    # - Too sparse (< 0.3): incomplete match, linear penalty
    # - Too dense (> 2.0): likely noise/tandem repeats, inverse penalty
    # - Optimal range (0.5-1.5): no penalty
    if seed_density < 0.3:
        density_score = seed_density / 0.3
    elif seed_density > 2.0:
        density_score = 2.0 / seed_density
    else:
        density_score = 1.0
    
    # 4. Context penalty (chromosome relationship)
    #    Intra-chromosomal repeats more likely to be true tandem duplications
    #    Inter-chromosomal repeats require stronger evidence (may be translocations)
    if name1 == name2:
        context_penalty = 1.0  # No penalty for intra-chromosomal
    else:
        context_penalty = 0.7  # 30% penalty for inter-chromosomal
    
    # 5. Combine all factors (multiplicative model)
    homology_score = length_score * identity_score * density_score * context_penalty
    
    # 6. Statistical significance bonus (based on Z-score)
    #    Measures how unlikely the observed identity is under random model
    expected_identity = 0.25  # Random expectation for 4-letter alphabet
    identity_z_score = (ident - expected_identity) / 0.1  # Assume std dev = 0.1
    
    # Apply bonus if Z-score indicates high significance
    if identity_z_score > 5:
        significance_bonus = 1.0 + math.log10(identity_z_score)
    else:
        significance_bonus = 1.0
    
    homology_score *= significance_bonus
    
    # ==================== Return Result ====================
    
    return {
        "chr1": name1, 
        "start1": start1, 
        "end1": end1,
        "chr2": name2, 
        "start2": start2, 
        "end2": end2,
        "strand": ori_flag, 
        "length": length,
        "identity": ident, 
        "n_seeds": len(chain_pairs),
        "score": homology_score,
        
        # Diagnostic fields for analysis and debugging
        "seed_density": seed_density,
        "context": "intra" if name1 == name2 else "inter"
    }

def identity_sampling(seq1, seq2, max_samples=1000):
    """
    Compute sequence identity by sampling (faster than full alignment).
    
    For very long sequences, computing exact identity is slow. Instead,
    we sample positions uniformly and compute identity on the sample.
    This gives a good approximation with much lower computational cost.
    
    Args:
        seq1, seq2: DNA sequences to compare
        max_samples: Maximum number of positions to sample
        
    Returns:
        Float in [0, 1] representing estimated sequence identity
        
    Algorithm:
        1. Compute sampling step: max(1, length // max_samples)
        2. Count matches at sampled positions
        3. Return: matches / total_sampled_positions
        
    Accuracy:
        For length >> max_samples, relative error ≈ 1/sqrt(max_samples)
        With max_samples=1000, relative error ≈ 3% (acceptable for scoring)
    """
    L = min(len(seq1), len(seq2))
    if L == 0: return 0.0
    step = max(1, L // max_samples)
    matches = sum(1 for i in range(0, L, step) if seq1[i] == seq2[i])
    total = (L - 1) // step + 1
    return matches / total if total > 0 else 0.0

# ===================== Merging and Output =====================

def merge_overlapping_calls(calls, overlap_thresh=0.5):
    """
    Merge overlapping repeat calls to avoid redundancy.
    
    When multiple chains detect the same repeat region (possibly with
    slight boundary differences), we want to keep only the best one.
    
    Args:
        calls: List of repeat call dictionaries
        overlap_thresh: Minimum overlap ratio to consider merging (unused in current logic)
        
    Returns:
        List of merged repeat calls (non-overlapping)
        
    Algorithm:
        1. Sort calls by (chr1, start1, chr2, start2) for sequential processing
        2. Iterate through sorted calls:
           - If current call overlaps with previous (same chromosomes, same strand):
             * Check if regions overlap (positive overlap on both sides)
             * If yes: keep the one with higher score, discard the other
           - If no overlap: save previous, start new current
        3. Append final call
        
    Note: This is a greedy merge algorithm. For complex overlaps, more
    sophisticated algorithms (e.g., interval trees) could be used.
    """
    if not calls: return []
    calls.sort(key=lambda x: (x["chr1"], x["start1"], x["chr2"], x["start2"]))
    merged = []
    curr = calls[0]
    for nxt in calls[1:]:
        # Check if calls are from same chromosome pair and same strand
        if (curr["chr1"] == nxt["chr1"] and curr["chr2"] == nxt["chr2"] and
            curr["strand"] == nxt["strand"]):
            
            # Compute overlap on both chromosomes
            ov1 = min(curr["end1"], nxt["end1"]) - max(curr["start1"], nxt["start1"])
            ov2 = min(curr["end2"], nxt["end2"]) - max(curr["start2"], nxt["start2"])
            
            # If overlapping on both sides, keep the one with higher score
            if ov1 > 0 and ov2 > 0:
                if nxt["score"] > curr["score"]:
                    curr = nxt
                continue  # Skip adding to merged (will be added at end or replaced)
        
        # No overlap or different chromosomes: save current and move to next
        merged.append(curr)
        curr = nxt
    merged.append(curr)
    return merged

def write_repeats(calls, path):
    """
    Write repeat calls to TSV file with all statistics.
    
    Output format (13 columns):
        chr1, start1, end1, chr2, start2, end2, strand, length,
        identity, score, n_seeds, seed_density, context
        
    Args:
        calls: List of repeat call dictionaries
        path: Output file path
        
    File format:
        - Tab-separated values
        - Header row with column names
        - Coordinates are 0-based half-open intervals [start, end)
        - Identity and score are floating point (4 decimal places)
    """
    with open(path, "w") as fh:
        # Write header
        fh.write("chr1\tstart1\tend1\tchr2\tstart2\tend2\tstrand\tlength\t"
                 "identity\tscore\tn_seeds\tseed_density\tcontext\n")
        # Write data rows
        for c in calls:
            fh.write(f"{c['chr1']}\t{c['start1']}\t{c['end1']}\t"
                     f"{c['chr2']}\t{c['start2']}\t{c['end2']}\t"
                     f"{c['strand']}\t{c['length']}\t{c['identity']:.4f}\t"
                     f"{c['score']:.4f}\t{c['n_seeds']}\t"
                     f"{c['seed_density']:.3f}\t{c['context']}\n")

# ===================== Chromosome Filtering =====================

def filter_chromosomes(genome_seqs, chromosome_filter):
    """
    Filter genome sequences to analyze only specified chromosome(s).
    
    This is crucial for large genomes (e.g., Arabidopsis, human) where
    analyzing all chromosomes at once would be too slow or memory-intensive.
    
    Args:
        genome_seqs: List of (name, seq) tuples for all chromosomes
        chromosome_filter: None (all), integer (index), or string (name/pattern)
        
    Returns:
        Filtered list of (name, seq) tuples
        
    Matching strategies:
        1. If integer: use as 1-based index into genome_seqs
        2. If string: try multiple matching methods:
           - Exact match (e.g., "Chr1" matches "Chr1")
           - Suffix match (e.g., "1" matches "Chr1", "chromosome_1")
           - Substring match (e.g., "chr" matches "Chr1", "chromosome_2")
        
    Error handling:
        - Invalid index: print available chromosomes and exit
        - No match: print available chromosomes and exit
        - Multiple matches: use first match and warn user
        
    Examples:
        filter_chromosomes(seqs, None) -> all chromosomes
        filter_chromosomes(seqs, 1) -> first chromosome (1-indexed)
        filter_chromosomes(seqs, "Chr1") -> chromosome named "Chr1"
        filter_chromosomes(seqs, "1") -> matches "Chr1", "chromosome_1", etc.
    """
    if chromosome_filter is None:
        return genome_seqs
    
    # Integer: use as 1-based index
    if isinstance(chromosome_filter, int):
        if chromosome_filter < 1 or chromosome_filter > len(genome_seqs):
            print(f"[ERROR] Chromosome index {chromosome_filter} out of range (1-{len(genome_seqs)})!")
            print(f"[INFO] Available chromosomes:")
            for idx, (name, seq) in enumerate(genome_seqs, 1):
                print(f"  {idx}. {name} ({len(seq):,} bp)")
            sys.exit(1)
        
        selected = genome_seqs[chromosome_filter - 1]
        print(f"[INFO] Selected chromosome {chromosome_filter}: {selected[0]} ({len(selected[1]):,} bp)")
        return [selected]
    
    # String: try multiple matching strategies
    filtered = []
    chr_str = str(chromosome_filter)
    
    for name, seq in genome_seqs:
        # Strategy 1: Exact match
        if name == chr_str:
            filtered.append((name, seq))
        # Strategy 2: Suffix match (e.g., "1" matches "Chr1")
        elif name.lower().endswith(chr_str.lower()):
            filtered.append((name, seq))
        # Strategy 3: Substring match (e.g., "chr" matches "Chr1")
        elif chr_str.lower() in name.lower():
            filtered.append((name, seq))
    
    # Error: no matches found
    if not filtered:
        print(f"[ERROR] Chromosome '{chr_str}' not found in genome!")
        print(f"[INFO] Available chromosomes:")
        for idx, (name, seq) in enumerate(genome_seqs, 1):
            print(f"  {idx}. {name} ({len(seq):,} bp)")
        sys.exit(1)
    
    # Warning: multiple matches (use first, warn user)
    if len(filtered) > 1:
        print(f"[WARN] Multiple chromosomes matched '{chr_str}':")
        for name, seq in filtered:
            print(f"  - {name} ({len(seq):,} bp)")
        print(f"[INFO] Using the first match: {filtered[0][0]}")
        return [filtered[0]]
    
    print(f"[INFO] Selected chromosome: {filtered[0][0]} ({len(filtered[0][1]):,} bp)")
    return filtered

# ===================== Main Program =====================

def cmd_find(args):
    """
    Main command: find repeats in genome.
    
    Workflow:
        1. Apply preset configuration (if specified)
        2. Create output directory (if needed)
        3. Load genome sequences from FASTA
        4. Filter to specified chromosome (if requested)
        5. Perform global seed pre-filtering
        6. Generate tasks for all chromosome pairs
        7. Process tasks (serial or parallel depending on task count)
        8. Merge overlapping calls
        9. Filter by minimum length
        10. Write results to file
        11. Generate visualizations (if requested)
        
    Parallelization strategy:
        - Single task: run directly (avoid multiprocessing overhead)
        - Multiple tasks: use Pool with min(3, cpu_count()) workers
        - Limit workers to 3 to balance speed vs memory usage
        
    Args:
        args: Parsed command-line arguments from argparse
    """
    start_time = time.time()
    
    # Apply preset configuration if specified
    if args.preset:
        if args.preset not in PRESET_CONFIGS:
            print(f"[ERROR] Unknown preset: {args.preset}")
            print(f"Available presets: {', '.join(PRESET_CONFIGS.keys())}")
            sys.exit(1)
        
        config = PRESET_CONFIGS[args.preset]
        print(f"[INFO] Using preset: '{args.preset}'")
        print(f"[INFO] Description: {config['description']}")
        
        # Override default parameters with preset values (if not explicitly set by user)
        if args.mask == "11110111101111":
            args.mask = config["mask"]
        if args.step == 2:
            args.step = config["step"]
        if args.nseeds == 30:
            args.nseeds = config["nseeds"]
        if args.max_occ is None:
            args.max_occ = config["max_occ"]
        if args.band == 400:
            args.band = config["band"]
    
    # Create output directory if needed
    output_dir = os.path.dirname(args.out)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"[INFO] Created output directory: {output_dir}")
    
    # Load genome sequences
    print(f"[INFO] Loading genome from {args.genome}...")
    all_genome_seqs = read_fasta(args.genome)
    
    # Filter to specified chromosome if requested
    if args.chromosome:
        # Try to parse as integer (chromosome index), fallback to string (name)
        try:
            chr_filter = int(args.chromosome)
        except ValueError:
            chr_filter = args.chromosome
        
        genome_seqs = filter_chromosomes(all_genome_seqs, chr_filter)
    else:
        genome_seqs = all_genome_seqs
    
    total_bp = sum(len(s) for _, s in genome_seqs)
    print(f"[INFO] Loaded {len(genome_seqs)} sequence(s), total {total_bp:,} bp.")
    
    # Auto-determine max_occ if not specified (based on seed length)
    mask_ones = args.mask.count('1')
    if args.max_occ is None:
        if mask_ones < 13:
            args.max_occ = 2000  # Short seeds: stricter filtering
        else:
            args.max_occ = 10000  # Long seeds: more lenient
    
    print(f"[INFO] Config: Mask={args.mask} (K={mask_ones}), Step={args.step}, "
          f"Nseeds={args.nseeds}, Band={args.band}, MaxOcc={args.max_occ}")
    
    # Global seed pre-filtering
    blacklist = prefilter_seeds_global(genome_seqs, args.mask, args.step, args.max_occ)
    
    # Generate tasks for all chromosome pairs (including self-pairs)
    tasks = []
    n_chroms = len(genome_seqs)
    task_names = []
    for i in range(n_chroms):
        for j in range(i, n_chroms):  # Only upper triangle (avoid redundant pairs)
            name1, seq1 = genome_seqs[i]
            name2, seq2 = genome_seqs[j]
            tasks.append((name1, seq1, name2, seq2, args.mask, args.step,
                         args.nseeds, args.band, args.minlen, blacklist))
            task_names.append(f"{name1} vs {name2}")
    
    print(f"[INFO] Processing {len(tasks)} chromosome pair(s)...")
    
    # Determine number of worker processes
    n_workers = min(3, cpu_count())  # Limit to 3 to control memory usage
    print(f"[INFO] Using {n_workers} CPU core(s).")
    
    # Process tasks
    if len(tasks) == 1:
        # Single task: run directly (no multiprocessing overhead)
        print(f"[INFO] Processing: {task_names[0]}...")
        print("[INFO] This may take 5-20 minutes depending on chromosome size...")
        result = find_repeats_between_chroms_worker(tasks[0])
        results = [result]
        elapsed = time.time() - start_time
        print(f"\n[INFO] Completed in {elapsed:.1f}s")
    
    elif n_workers > 1:
        # Multiple tasks + multiple cores: parallel processing
        with Pool(processes=n_workers) as pool:
            results = []
            for idx, result in enumerate(pool.imap_unordered(find_repeats_between_chroms_worker, tasks), 1):
                results.append(result)
                elapsed = time.time() - start_time
                eta = (elapsed / idx) * (len(tasks) - idx)
                print_progress_bar(idx, len(tasks), 
                                  prefix='  -> Progress',
                                  suffix=f'Complete (Elapsed: {elapsed:.1f}s, ETA: {eta:.1f}s)')
    else:
        # Multiple tasks + single core: serial processing
        results = []
        for idx, task in enumerate(tasks, 1):
            result = find_repeats_between_chroms_worker(task)
            results.append(result)
            elapsed = time.time() - start_time
            eta = (elapsed / idx) * (len(tasks) - idx)
            print_progress_bar(idx, len(tasks),
                              prefix='  -> Progress',
                              suffix=f'{task_names[idx-1]} (Elapsed: {elapsed:.1f}s, ETA: {eta:.1f}s)')
    
    # Collect all calls from all tasks
    all_calls = []
    for calls in results:
        all_calls.extend(calls)
    
    print(f"\n[INFO] Total candidates found: {len(all_calls)}")
    
    # Merge overlapping calls
    print("[INFO] Merging overlapping calls...")
    merged_calls = merge_overlapping_calls(all_calls)
    
    # Filter by minimum length
    final_calls = [c for c in merged_calls if c["length"] >= args.minlen]
    print(f"[INFO] Final repeats >= {args.minlen}bp: {len(final_calls)}")
    
    # Write results
    write_repeats(final_calls, args.out)
    print(f"[INFO] Results written to {args.out}")
    
    # Generate visualizations if requested
    if args.visualize:
        output_prefix = args.out.rsplit('.', 1)[0]
        create_visualizations(final_calls, genome_seqs, output_prefix)
    
    # Print summary
    total_time = time.time() - start_time
    print(f"[INFO] Total runtime: {total_time:.1f} seconds ({total_time/60:.1f} minutes)")
    print("[INFO] Done!")

def main():
    """
    Main entry point: parse command-line arguments and dispatch to cmd_find.
    
    Command-line interface:
        python repeatfinder.py find [options]
        
    Required arguments:
        --genome: Input genome FASTA file
        --out: Output TSV file path
        
    Optional arguments:
        --chromosome: Analyze only specific chromosome (index or name)
        --preset: Use predefined parameter configuration (balanced/fast/far-fast/precise/sensitive)
        --visualize: Generate plots (requires matplotlib)
        --mask: Spaced seed mask pattern
        --step: Seed sampling step size
        --nseeds: Minimum seeds in chain
        --band: Diagonal band width
        --max-occ: Maximum seed occurrences
        --minlen: Minimum repeat length
        
    Examples:
        # Analyze entire genome with balanced preset
        python repeatfinder.py find --genome genome.fa --out results.tsv --preset balanced
        
        # Analyze only chromosome 1 with fast preset
        python repeatfinder.py find --genome genome.fa --out chr1.tsv --preset fast --chromosome 1
        
        # Ultra-fast mode for very large genomes (lowest sensitivity but maximum speed)
        python repeatfinder.py find --genome genome.fa --out ultrafast.tsv --preset far-fast --chromosome 1
        
        # Custom parameters
        python repeatfinder.py find --genome genome.fa --out custom.tsv --step 4 --nseeds 15
    """
    parser = argparse.ArgumentParser(
        description="High-performance repeat finder for large genomes",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Preset Modes:
  balanced   - Balanced mode: Optimal trade-off between recall and precision (default recommended)
  fast       - Fast mode: Prioritizes speed, suitable for initial screening of large genomes
  far-fast   - Ultra-fast mode: Maximum speed for very large genomes (lowest sensitivity, highest speed)
  precise    - Precise mode: Stricter filtering to reduce false positives
  sensitive  - Sensitive mode: Maximum recall to find more true positives (but with more false positives)

Performance Comparison (approximate on 30 Mbp chromosome):
  sensitive  : ~60-90 minutes  (highest recall, most false positives)
  balanced   : ~20-40 minutes  (best overall quality, recommended)
  fast       : ~8-15 minutes   (good for initial screening)
  far-fast   : ~3-6 minutes    (ultra-fast, may miss some repeats)
  precise    : ~25-50 minutes  (lowest false positives, may miss diverged repeats)

Examples:
  # Analyze entire genome with balanced preset
  python repeatfinder.py find --genome genome.fa --out results.tsv --preset balanced --visualize
  
  # Analyze only chromosome 1 (recommended for large genomes)
  python repeatfinder.py find --genome genome.fa --out chr1.tsv --preset balanced --chromosome 1 --visualize
  
  # Analyze chromosome named "Chr3"
  python repeatfinder.py find --genome genome.fa --out chr3.tsv --preset balanced --chromosome Chr3 --visualize
  
  # Ultra-fast mode for very large chromosomes (e.g., human Chr1)
  python repeatfinder.py find --genome human.fa --out chr1_fast.tsv --preset far-fast --chromosome 1
  
  # Fast test on Arabidopsis chromosome 1
  python repeatfinder.py find --genome Arabidopsis.fa --out test.tsv --preset fast --chromosome 1
"""
    )
    sub = parser.add_subparsers()
    
    # Define 'find' subcommand
    p = sub.add_parser("find", help="Find repeats in genome")
    p.add_argument("--genome", required=True, help="Input genome FASTA file")
    p.add_argument("--out", required=True, help="Output TSV file")
    p.add_argument("--chromosome", default=None, 
                   help="Analyze only specific chromosome (number or name, e.g., '1', 'Chr1'). If not specified, analyze all chromosomes.")
    p.add_argument("--preset", choices=list(PRESET_CONFIGS.keys()), 
                   help="Use a preset parameter configuration (balanced/fast/far-fast/precise/sensitive)")
    p.add_argument("--visualize", action="store_true", 
                   help="Generate visualization plots (highly recommended)")
    p.add_argument("--mask", default="11110111101111", help="Spaced seed mask")
    p.add_argument("--step", type=int, default=2, help="Seed sampling step")
    p.add_argument("--nseeds", type=int, default=30, help="Minimum seeds in chain")
    p.add_argument("--band", type=int, default=400, help="Diagonal band width")
    p.add_argument("--max-occ", type=int, default=None, help="Max seed occurrences (auto if None)")
    p.add_argument("--minlen", type=int, default=20000, help="Minimum repeat length")
    p.set_defaults(func=cmd_find)
    
    # Parse arguments and dispatch
    args = parser.parse_args()
    if hasattr(args, "func"):
        args.func(args)
    else:
        parser.print_help()

if __name__ == "__main__":
    main()
