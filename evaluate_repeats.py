#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
evaluate_repeats.py (Updated to support minimap2 format)

Compare ground_truth_repeats.tsv, minimap2 output, and repeatfinder.py output
to evaluate the accuracy of repeat detection algorithms.

This script supports multiple evaluation scenarios:
1. Your algorithm vs ground truth (gold standard evaluation)
2. Your algorithm vs minimap2 (consistency check against established tool)
3. Minimap2 vs ground truth (baseline performance measurement)

Usage examples:
  # Scenario 1: Evaluate your results against ground truth
  python evaluate_repeats.py \
      --truth ground_truth_repeats.tsv \
      --pred syn_my_result.tsv \
      --overlap 0.8
  
  # Scenario 2: Compare your results with minimap2 (treat minimap2 as reference)
  python evaluate_repeats.py \
      --truth syn_adjust.tsv \
      --pred syn_my_result.tsv \
      --overlap 0.5 \
      --truth-format minimap
  
  # Scenario 3: Evaluate minimap2 against ground truth (baseline benchmark)
  python evaluate_repeats.py \
      --truth ground_truth_repeats.tsv \
      --pred syn_adjust.tsv \
      --overlap 0.8 \
      --pred-format minimap
"""

import argparse

def read_truth_repeats(path):
    """
    Read ground truth repeat file in standard format.
    
    Expected file format (TSV with header):
      repeat_id  source_start_1based  source_end_1based  target_start_1based  target_end_1based  strand  ...
    
    This format is produced by data.py's build_dataset() function and represents
    the "true" repeat locations that were artificially inserted into the synthetic genome.
    
    Args:
        path: Path to ground truth TSV file
        
    Returns:
        List of tuples [(chr1, start1, end1, chr2, start2, end2), ...]
        where coordinates are in 0-based half-open format [start, end)
        
    Coordinate conversion:
        - Input file uses 1-based inclusive intervals (biological convention)
        - Output uses 0-based half-open intervals (Python/programming convention)
        - Conversion: start_0based = start_1based - 1, end_0based = end_1based (no change)
        
    Example:
        Input line: "1  1000  2000  5000  6000  +  ..."
        Output: ("chr1", 999, 2000, "chr1", 4999, 6000)
        
    Assumptions:
        - All repeats are on chromosome "chr1" (single-chromosome test dataset)
        - First line is header (skipped)
        - Minimum 5 columns (repeat_id, source_start, source_end, target_start, target_end)
    """
    reps = []
    with open(path, "r") as fh:
        header_skipped = False
        for line in fh:
            line = line.strip()
            # Skip empty lines and comment lines
            if not line or line.startswith("#"):
                continue
            
            cols = line.split("\t")
            
            # Skip header row (first non-comment line)
            if not header_skipped:
                header_skipped = True
                print(f"[DEBUG] Skipping header: {cols[0]}")
                continue
            
            # Validate minimum column count
            if len(cols) < 5:
                continue
            
            try:
                # Parse columns (0-indexed):
                # 0 = repeat_id (string identifier)
                # 1 = source_start (1-based, inclusive)
                # 2 = source_end (1-based, inclusive)
                # 3 = target_start (1-based, inclusive)
                # 4 = target_end (1-based, inclusive)
                # 5 = strand (optional: "+" or "-")
                repeat_id = cols[0]
                
                # Convert from 1-based inclusive to 0-based half-open
                # Example: [1000, 2000] (1-based inclusive) -> [999, 2000) (0-based half-open)
                source_start = int(cols[1]) - 1  # Subtract 1 for 0-based start
                source_end = int(cols[2])        # No change: 1-based inclusive end = 0-based exclusive end
                target_start = int(cols[3]) - 1  # Subtract 1 for 0-based start
                target_end = int(cols[4])        # No change
                
                strand = cols[5] if len(cols) > 5 else "+"
                
                # Assume all repeats are on chr1 (single-chromosome synthetic genome)
                reps.append(("chr1", source_start, source_end, "chr1", target_start, target_end))
                
                print(f"[DEBUG] Truth repeat #{repeat_id}: "
                      f"source=[{source_start}, {source_end}), len={source_end - source_start}, "
                      f"target=[{target_start}, {target_end}), len={target_end - target_start}, "
                      f"strand={strand}")
                
            except (ValueError, IndexError) as e:
                print(f"[WARN] Failed to parse line: {line[:50]}... Error: {e}")
                continue
    
    return reps


def read_minimap_repeats(path, minlen=15000):
    """
    Read minimap2 self-alignment output in PAF (Pairwise mApping Format).
    
    PAF is the standard output format of minimap2 for sequence alignments.
    When aligning a genome to itself, high-quality alignments (excluding the
    trivial self-match) represent repeat regions.
    
    Expected file format (TSV with header):
      qname  qlen  qstart  qend  strand  tname  tlen  tstart  tend  n_match  aln_len  mapq
    
    Args:
        path: Path to minimap2 PAF file
        minlen: Minimum repeat length to keep (default: 15000 bp)
               Repeats shorter than this are filtered out to reduce noise
               
    Returns:
        List of tuples [(chr1, start1, end1, chr2, start2, end2), ...]
        where coordinates are in 0-based half-open format [start, end)
        
    Coordinate system:
        - PAF format already uses 0-based half-open intervals [start, end)
        - No conversion needed (unlike ground truth format)
        
    Filtering criteria:
        1. Self-matches: Skip alignments where qstart ≈ tstart (within 100bp)
           These are trivial matches of a region to itself
           
        2. Length filter: Skip alignments shorter than minlen
           Short repeats are more likely to be random similarities
           
        3. Header: Skip first line (column names)
        
    Example:
        Input line: "chr1  3000000  10000  30000  +  chr1  3000000  50000  70000  19500  20000  60"
        Output: ("chr1", 10000, 30000, "chr1", 50000, 70000)
        Interpretation: 20kb repeat between positions 10-30kb and 50-70kb
        
    Quality metrics (for reference, not used in filtering):
        - n_match: Number of matching bases in alignment
        - aln_len: Total alignment length (matches + mismatches + gaps)
        - mapq: Mapping quality (Phred-scaled probability of incorrect alignment)
    """
    reps = []
    with open(path, "r") as fh:
        header_skipped = False
        for line_num, line in enumerate(fh, start=1):
            line = line.strip()
            # Skip empty lines and comment lines
            if not line or line.startswith("#"):
                continue
            
            cols = line.split("\t")
            
            # Skip header row
            if not header_skipped:
                header_skipped = True
                print(f"[DEBUG] Skipping minimap header: {cols[0]}")
                continue
            
            # Validate minimum column count (need at least 9 columns)
            if len(cols) < 9:
                continue
            
            try:
                # Parse PAF columns:
                qname = cols[0]    # Query sequence name
                qstart = int(cols[2])  # Query start (0-based)
                qend = int(cols[3])    # Query end (0-based, exclusive)
                strand = cols[4]       # Strand: "+" or "-"
                tname = cols[5]        # Target sequence name
                tstart = int(cols[7])  # Target start (0-based)
                tend = int(cols[8])    # Target end (0-based, exclusive)
                
                # Filter 1: Skip trivial self-matches
                # A region should not match itself at the same position
                # Allow small tolerance (100bp) for alignment boundary variations
                if qname == tname and abs(qstart - tstart) < 100 and abs(qend - tend) < 100:
                    continue
                
                # Calculate repeat lengths
                q_length = qend - qstart
                t_length = tend - tstart
                length = min(q_length, t_length)  # Use shorter length (conservative)
                
                # Filter 2: Skip short repeats (likely random similarity)
                if length < minlen:
                    continue
                
                # Add to results
                reps.append((qname, qstart, qend, tname, tstart, tend))
                
                print(f"[DEBUG] Minimap repeat #{line_num}: "
                      f"{qname}:[{qstart}, {qend}) -> {tname}:[{tstart}, {tend}), "
                      f"len={length}")
                
            except (ValueError, IndexError) as e:
                print(f"[WARN] Line {line_num}: Failed to parse - {e}")
                continue
    
    return reps


def read_pred_repeats(path):
    """
    Read predicted repeat file from repeatfinder.py.
    
    Expected file format (TSV with header):
      chr1  start1  end1  chr2  start2  end2  strand  length  identity  score  n_seeds  seed_density  context
    
    This is the output format of repeatfinder.py's write_repeats() function.
    
    Args:
        path: Path to predicted repeats TSV file
        
    Returns:
        List of tuples [(chr1, start1, end1, chr2, start2, end2), ...]
        where coordinates are in 0-based half-open format [start, end)
        
    Coordinate system:
        - repeatfinder.py already outputs 0-based half-open intervals [start, end)
        - No conversion needed
        
    Note: Only first 6 columns are used for evaluation (chromosome names and coordinates)
          Other columns (strand, length, identity, score, etc.) are for analysis but not
          required for overlap calculation.
          
    Example:
        Input line: "chr1  10000  30000  chr1  50000  70000  +  20000  0.95  85.4  45  0.85  intra"
        Output: ("chr1", 10000, 30000, "chr1", 50000, 70000)
        Interpretation: Predicted 20kb repeat between [10000, 30000) and [50000, 70000)
    """
    reps = []
    with open(path, "r") as fh:
        header_skipped = False
        for line in fh:
            line = line.strip()
            # Skip empty lines and comment lines
            if not line or line.startswith("#"):
                continue
            
            cols = line.split("\t")
            
            # Skip header row
            if not header_skipped:
                header_skipped = True
                print(f"[DEBUG] Skipping pred header: {cols[0]}")
                continue
            
            # Validate minimum column count (need at least 6: chr1, start1, end1, chr2, start2, end2)
            if len(cols) < 6:
                continue
            
            try:
                # Parse coordinates (already 0-based half-open from repeatfinder.py)
                chr1 = cols[0]
                start1 = int(cols[1])
                end1 = int(cols[2])
                chr2 = cols[3]
                start2 = int(cols[4])
                end2 = int(cols[5])
                
                # Add to results
                reps.append((chr1, start1, end1, chr2, start2, end2))
                
                print(f"[DEBUG] Pred repeat: "
                      f"{chr1}:[{start1}, {end1}) -> {chr2}:[{start2}, {end2}), "
                      f"len1={end1 - start1}, len2={end2 - start2}")
                
            except (ValueError, IndexError) as e:
                print(f"[WARN] Failed to parse line: {line[:50]}... Error: {e}")
                continue
    
    return reps


def interval_overlap_ratio(a_start, a_end, b_start, b_end):
    """
    Calculate overlap ratio between two intervals using Jaccard-like metric.
    
    This function measures how well two intervals overlap by computing:
        overlap_ratio = intersection_length / min(length_a, length_b)
        
    Rationale for using min() instead of union (Jaccard):
        - If interval A is [0, 100) and B is [50, 150), intersection = 50bp
        - Jaccard would give 50/150 = 0.33 (low score)
        - Min-based gives 50/100 = 0.5 (more reasonable)
        - This is more forgiving when one prediction is slightly longer than the truth
        
    Args:
        a_start, a_end: First interval [a_start, a_end) (0-based, half-open)
        b_start, b_end: Second interval [b_start, b_end) (0-based, half-open)
        
    Returns:
        Float in [0, 1] representing overlap ratio
        - 0.0: No overlap
        - 1.0: Perfect overlap (one interval completely contains the other)
        - Values between 0 and 1: Partial overlap
        
    Algorithm:
        1. Compute intersection: max(0, min(a_end, b_end) - max(a_start, b_start))
        2. Compute lengths: a_length = a_end - a_start, b_length = b_end - b_start
        3. Return: intersection / min(a_length, b_length)
        
    Edge cases:
        - If either interval has length 0, return 0.0
        - If intervals don't overlap, intersection = 0, return 0.0
        
    Example:
        interval_overlap_ratio(0, 100, 50, 150) = 50 / min(100, 100) = 0.5
        interval_overlap_ratio(0, 100, 200, 300) = 0 / 100 = 0.0 (no overlap)
        interval_overlap_ratio(0, 100, 25, 75) = 50 / min(100, 50) = 1.0 (complete containment)
    """
    # Compute intersection length (max with 0 to handle non-overlapping intervals)
    inter = max(0, min(a_end, b_end) - max(a_start, b_start))
    
    # Compute interval lengths (max with 0 to handle invalid intervals)
    len_a = max(0, a_end - a_start)
    len_b = max(0, b_end - b_start)
    
    # Handle edge case: if either interval has zero length, no meaningful overlap
    if len_a == 0 or len_b == 0:
        return 0.0
    
    # Return ratio: intersection / smaller interval length
    return inter / min(len_a, len_b)


def evaluate_truth_vs_pred(truth, pred, overlap_thresh=0.8):
    """
    Evaluate predicted repeats against ground truth using greedy matching.
    
    This function implements a greedy matching algorithm to pair each ground truth
    repeat with at most one predicted repeat, maximizing overall overlap quality.
    
    Matching criteria:
        1. Chromosomes must match (chr1 == chr1, chr2 == chr2)
        2. Both intervals must have overlap_ratio >= overlap_thresh
        3. Use minimum of two overlap ratios as overall score
        4. Each prediction can only match one truth (prevent double-counting)
        
    Args:
        truth: List of ground truth repeats [(chr1, s1, e1, chr2, s2, e2), ...]
        pred: List of predicted repeats [(chr1, s1, e1, chr2, s2, e2), ...]
        overlap_thresh: Minimum overlap ratio required for match (default: 0.8 = 80%)
                       Lower values are more lenient, higher values more strict
                       
    Returns:
        Tuple of (tp, fp, fn, precision, recall, f1, matched_details):
        - tp: True positives (# truth repeats successfully matched)
        - fp: False positives (# predictions not matching any truth)
        - fn: False negatives (# truth repeats not matched by any prediction)
        - precision: tp / (tp + fp), measures exactness
        - recall: tp / (tp + fn), measures completeness
        - f1: Harmonic mean of precision and recall (2*p*r / (p+r))
        - matched_details: List of (truth_idx, pred_idx, truth, pred, score) tuples
        
    Algorithm:
        For each ground truth repeat:
            1. Search all unmatched predictions
            2. Compute overlap score for matching chromosomes
            3. Keep track of best-scoring prediction (if score >= threshold)
            4. If found, mark as matched (prevent reuse)
            5. If not found, count as false negative
        
        After processing all truth repeats:
            - tp = number of matched truth repeats
            - fp = number of unmatched predictions
            - fn = number of unmatched truth repeats
            
    Matching strategy (greedy):
        - Process truth repeats in order
        - For each truth, pick best available prediction
        - This is not globally optimal but fast and simple
        - More sophisticated algorithms (e.g., Hungarian algorithm) could
          find globally optimal matching but are unnecessary for this use case
          
    Example scenario:
        Truth: 20 repeats
        Predictions: 18 repeats
        Matches found: 15 (with overlap >= 0.8)
        
        Results:
        - TP = 15 (15 truth repeats successfully matched)
        - FP = 3 (18 - 15 unmatched predictions, likely false positives)
        - FN = 5 (20 - 15 unmatched truth repeats, algorithm missed these)
        - Precision = 15/18 = 0.833 (83.3% of predictions are correct)
        - Recall = 15/20 = 0.75 (found 75% of true repeats)
        - F1 = 2*0.833*0.75 / (0.833+0.75) = 0.789
    """
    matched_pred_indices = set()  # Track which predictions have been used
    matched_details = []          # Store details of successful matches
    tp = 0                        # True positives counter

    print(f"\n[INFO] Matching with overlap threshold = {overlap_thresh}")
    print("=" * 80)

    # Process each ground truth repeat
    for i, t in enumerate(truth):
        t_chr1, t_s1, t_e1, t_chr2, t_s2, t_e2 = t
        best_j = None          # Index of best matching prediction
        best_score = 0.0       # Score of best match
        
        print(f"\n[Truth #{i+1}] {t_chr1}:[{t_s1}, {t_e1}) -> {t_chr2}:[{t_s2}, {t_e2})")
        print(f"  Searching for matches in {len(pred)} predictions...")
        
        # Search through all predictions for best match
        for j, p in enumerate(pred):
            # Skip predictions that have already been matched
            if j in matched_pred_indices:
                continue
                
            p_chr1, p_s1, p_e1, p_chr2, p_s2, p_e2 = p
            
            # Chromosome names must match exactly
            if t_chr1 != p_chr1 or t_chr2 != p_chr2:
                continue
            
            # Compute overlap ratios for both intervals
            # r1: overlap ratio for first interval (source region)
            # r2: overlap ratio for second interval (target region)
            r1 = interval_overlap_ratio(t_s1, t_e1, p_s1, p_e1)
            r2 = interval_overlap_ratio(t_s2, t_e2, p_s2, p_e2)
            
            # Overall score is minimum of two ratios (both must be good)
            # This ensures that BOTH ends of the repeat are well-predicted
            score = min(r1, r2)
            
            # Print candidates that meet threshold (for debugging)
            if score >= overlap_thresh:
                print(f"    Candidate pred #{j}: overlap_r1={r1:.3f}, overlap_r2={r2:.3f}, score={score:.3f}")
            
            # Update best match if this is better and meets threshold
            if score >= overlap_thresh and score > best_score:
                best_score = score
                best_j = j
        
        # Process matching result for this truth repeat
        if best_j is not None:
            # Match found: increment TP, mark prediction as used
            tp += 1
            matched_pred_indices.add(best_j)
            p = pred[best_j]
            matched_details.append((i, best_j, t, p, best_score))
            print(f"  ✓ MATCHED with pred #{best_j} (score={best_score:.3f})")
            print(f"    Pred: {p[0]}:[{p[1]}, {p[2]}) -> {p[3]}:[{p[4]}, {p[5]})")
        else:
            # No match found: this truth repeat will count as FN
            print(f"  ✗ NOT MATCHED (no prediction meets threshold)")

    # Calculate final statistics
    # False positives: predictions that didn't match any truth
    fp = len(pred) - tp
    
    # False negatives: truth repeats that weren't matched by any prediction
    fn = len(truth) - tp

    # Precision: what fraction of predictions are correct?
    prec = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    
    # Recall: what fraction of truth repeats were found?
    rec  = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    
    # F1 score: harmonic mean of precision and recall
    # Balances both metrics (neither can be ignored)
    f1   = 2 * prec * rec / (prec + rec) if (prec + rec) > 0 else 0.0

    return tp, fp, fn, prec, rec, f1, matched_details


def main():
    """
    Main entry point: parse arguments and run evaluation.
    
    This script supports flexible evaluation scenarios:
    1. Different input formats (ground truth, repeatfinder, minimap2)
    2. Adjustable overlap thresholds (stricter or more lenient matching)
    3. Length filtering for minimap2 input (reduce noise)
    
    Typical workflows:
    
    A) Validate your algorithm against synthetic ground truth:
       python evaluate_repeats.py \
           --truth ground_truth_repeats.tsv \
           --pred my_results.tsv \
           --overlap 0.8
       
       Interpretation: How well does your algorithm find known repeats?
       
    B) Compare your algorithm with minimap2:
       python evaluate_repeats.py \
           --truth minimap2_results.tsv \
           --pred my_results.tsv \
           --overlap 0.5 \
           --truth-format minimap
       
       Interpretation: How similar are your results to a mature tool?
       
    C) Benchmark minimap2 against ground truth:
       python evaluate_repeats.py \
           --truth ground_truth_repeats.tsv \
           --pred minimap2_results.tsv \
           --overlap 0.8 \
           --pred-format minimap
       
       Interpretation: What's the baseline performance to beat?
       
    Output interpretation:
        - High precision, low recall: algorithm is conservative (few false positives but misses many)
        - Low precision, high recall: algorithm is aggressive (finds most but many false positives)
        - High precision and recall: algorithm works well
        - F1 score: single metric balancing both (higher is better, max = 1.0)
    """
    ap = argparse.ArgumentParser(
        description="Evaluate long repeat predictions against ground truth or minimap2."
    )
    ap.add_argument("--truth", required=True, 
                    help="Ground truth or reference repeats TSV file. "
                         "Format depends on --truth-format.")
    ap.add_argument("--pred", required=True, 
                    help="Predicted repeats TSV file to evaluate. "
                         "Format depends on --pred-format.")
    ap.add_argument("--overlap", type=float, default=0.8,
                    help="Overlap threshold for matching (default: 0.8 = 80%%). "
                         "Lower values are more lenient, higher values more strict. "
                         "Typical range: 0.5 (loose) to 0.9 (very strict).")
    ap.add_argument("--truth-format", choices=["groundtruth", "minimap"], default="groundtruth",
                    help="Format of truth file: "
                         "'groundtruth' (default, from data.py) or "
                         "'minimap' (PAF format from minimap2).")
    ap.add_argument("--pred-format", choices=["repeatfinder", "minimap"], default="repeatfinder",
                    help="Format of pred file: "
                         "'repeatfinder' (default, from repeatfinder.py) or "
                         "'minimap' (PAF format from minimap2).")
    ap.add_argument("--minlen", type=int, default=15000,
                    help="Minimum repeat length for minimap format (default: 15000 bp). "
                         "Shorter repeats are filtered out to reduce noise.")
    args = ap.parse_args()

    # Step 1: Load reference/truth repeats
    print("=" * 80)
    print("STEP 1: Loading truth/reference repeats")
    print("=" * 80)
    
    if args.truth_format == "minimap":
        truth = read_minimap_repeats(args.truth, minlen=args.minlen)
    else:
        truth = read_truth_repeats(args.truth)
    
    print(f"\n[INFO] Loaded {len(truth)} reference repeats from {args.truth}")

    # Step 2: Load predicted repeats
    print("\n" + "=" * 80)
    print("STEP 2: Loading predicted repeats")
    print("=" * 80)
    
    if args.pred_format == "minimap":
        pred = read_minimap_repeats(args.pred, minlen=args.minlen)
    else:
        pred = read_pred_repeats(args.pred)
    
    print(f"\n[INFO] Loaded {len(pred)} predicted repeats from {args.pred}")

    # Validate that data was loaded successfully
    if len(truth) == 0:
        print("\n[ERROR] No reference repeats loaded! Check file format.")
        return
    
    if len(pred) == 0:
        print("\n[ERROR] No predicted repeats loaded! Check file format.")
        return

    # Step 3: Evaluate predictions
    print("\n" + "=" * 80)
    print("STEP 3: Evaluating predictions")
    print("=" * 80)
    
    tp, fp, fn, prec, rec, f1, matches = evaluate_truth_vs_pred(
        truth, pred, overlap_thresh=args.overlap
    )

    # Step 4: Print final results
    print("\n" + "=" * 80)
    print("FINAL RESULTS")
    print("=" * 80)
    print(f"TP (True Positives)  = {tp}")
    print(f"FP (False Positives) = {fp}")
    print(f"FN (False Negatives) = {fn}")
    print(f"Precision = {prec:.3f}")
    print(f"Recall    = {rec:.3f}")
    print(f"F1        = {f1:.3f}")
    print("=" * 80)

    # Print detailed match information (optional, for debugging)
    if matches:
        print("\nMatched pairs:")
        for truth_idx, pred_idx, t, p, score in matches:
            print(f"  Truth #{truth_idx+1} <-> Pred #{pred_idx} (overlap score={score:.3f})")


if __name__ == "__main__":
    main()
