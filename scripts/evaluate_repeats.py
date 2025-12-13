import argparse

###############################################################################
#                               PARSE TRUTH FILE
###############################################################################

def read_truth_repeats(path):
    """
    Read synthetic dataset ground truth repeats.
    Expected columns:
    repeat_id  source_start_1based  source_end_1based  target_start_1based  target_end_1based  strand ...
    """
    reps = []
    with open(path, "r") as fh:
        header_skipped = False
        
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            
            cols = line.split("\t")

            # Skip header line
            if not header_skipped:
                header_skipped = True
                print(f"[DEBUG] Skipping ground truth header: {cols[0]}")
                continue
            
            if len(cols) < 5:
                continue

            try:
                repeat_id = cols[0]

                # Convert 1-based inclusive → 0-based half-open
                src_start = int(cols[1]) - 1
                src_end   = int(cols[2])
                tgt_start = int(cols[3]) - 1
                tgt_end   = int(cols[4])

                reps.append(("chr1", src_start, src_end, "chr1", tgt_start, tgt_end))

                print(f"[DEBUG] Truth #{repeat_id}: "
                      f"src=[{src_start},{src_end}), tgt=[{tgt_start},{tgt_end})")

            except Exception as e:
                print(f"[WARN] Failed to parse truth: {line}  Error={e}")

    return reps


###############################################################################
#                         PARSE MINIMAP2 (PAF) OUTPUT
###############################################################################

def read_minimap_repeats(path, minlen=15000):
    """
    Parse minimap2 PAF output (only first 9 columns used).
    """
    reps = []
    with open(path, "r") as fh:
        header_skipped = False
        for line_num, line in enumerate(fh, start=1):
            line = line.strip()
            if not line:
                continue

            cols = line.split("\t")

            # Skip header
            if not header_skipped:
                header_skipped = True
                print(f"[DEBUG] Skipping minimap header: {cols[0]}")
                continue

            if len(cols) < 9:
                continue

            try:
                qname  = cols[0]
                qstart = int(cols[2])
                qend   = int(cols[3])
                strand = cols[4]
                tname  = cols[5]
                tstart = int(cols[7])
                tend   = int(cols[8])

                # Remove trivial self-match
                if qname == tname and abs(qstart - tstart) < 100 and abs(qend - tend) < 100:
                    continue

                length = min(qend - qstart, tend - tstart)
                if length < minlen:
                    continue

                reps.append((qname, qstart, qend, tname, tstart, tend))

                print(f"[DEBUG] Minimap repeat #{line_num}: "
                      f"{qname}:[{qstart},{qend}) -> {tname}:[{tstart},{tend}), len={length}")

            except Exception as e:
                print(f"[WARN] Minimpa2 parse error line {line_num}: {e}")

    return reps


###############################################################################
#                          PARSE MUMMER show-coords OUTPUT
###############################################################################

def read_mummer_repeats(path, minlen=15000):
    """
    Parse MUMmer show-coords -rcl TSV.
    Coordinates are 1-based inclusive.
    """
    reps = []
    header_skipped = False

    with open(path, "r") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue

            # Skip NUCMER header lines
            if line.startswith("/") or line.startswith("NUCMER"):
                continue

            cols = line.split("\t")

            # Skip header line starting with [S1]
            if not header_skipped:
                if cols[0] == "[S1]":
                    header_skipped = True
                continue

            if len(cols) < 12:
                continue

            try:
                S1 = int(cols[0])
                E1 = int(cols[1])
                S2 = int(cols[2])
                E2 = int(cols[3])
                LEN1 = int(cols[4])
                chrom1 = cols[11]
                chrom2 = cols[12] if len(cols) > 12 else chrom1

                if LEN1 < minlen:
                    continue

                # Convert to 0-based half-open
                reps.append((chrom1, S1 - 1, E1, chrom2, S2 - 1, E2))

                print(f"[DEBUG] MUMmer repeat: {chrom1}:[{S1-1},{E1}) -> {chrom2}:[{S2-1},{E2}), len={LEN1}")

            except Exception:
                continue

    return reps


###############################################################################
#                               PARSE REPEATFINDER OUTPUT
###############################################################################

def read_pred_repeats(path):
    """
    Parse repeatfinder output: 6 columns (chr1, start1, end1, chr2, start2, end2)
    Already in 0-based half-open.
    """
    reps = []
    header_skipped = False

    with open(path, "r") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue

            cols = line.split("\t")

            # Skip header
            if not header_skipped:
                header_skipped = True
                print(f"[DEBUG] Skipping pred header: {cols[0]}")
                continue

            if len(cols) < 6:
                continue

            try:
                chr1, s1, e1, chr2, s2, e2 = cols[:6]
                reps.append((chr1, int(s1), int(e1), chr2, int(s2), int(e2)))
            except Exception:
                continue

    return reps


###############################################################################
#                                OVERLAP LOGIC
###############################################################################

def interval_overlap_ratio(a_start, a_end, b_start, b_end):
    inter = max(0, min(a_end, b_end) - max(a_start, b_start))
    len_a = max(0, a_end - a_start)
    len_b = max(0, b_end - b_start)
    if len_a == 0 or len_b == 0:
        return 0.0
    return inter / min(len_a, len_b)


###############################################################################
#                          MATCHING TRUTH VS PREDICTED
###############################################################################

def evaluate_truth_vs_pred(truth, pred, overlap_thresh=0.8):
    matched_preds = set()
    matched_details = []

    tp = 0

    print(f"\n[INFO] Matching (overlap threshold = {overlap_thresh})")
    print("=" * 80)

    for i, t in enumerate(truth):
        t_chr1, t_s1, t_e1, t_chr2, t_s2, t_e2 = t
        best_j = None
        best_score = 0.0

        print(f"\n[Truth #{i+1}] {t_chr1}:[{t_s1},{t_e1}) -> {t_chr2}:[{t_s2},{t_e2})")

        for j, p in enumerate(pred):
            if j in matched_preds:
                continue

            p_chr1, p_s1, p_e1, p_chr2, p_s2, p_e2 = p

            if p_chr1 != t_chr1 or p_chr2 != t_chr2:
                continue

            r1 = interval_overlap_ratio(t_s1, t_e1, p_s1, p_e1)
            r2 = interval_overlap_ratio(t_s2, t_e2, p_s2, p_e2)
            score = min(r1, r2)

            if score >= overlap_thresh and score > best_score:
                best_score = score
                best_j = j

        if best_j is not None:
            tp += 1
            matched_preds.add(best_j)
            matched_details.append((i, best_j, t, pred[best_j], best_score))
            print(f"  ✓ MATCHED with pred #{best_j}, score={best_score:.3f}")
        else:
            print(f"  ✗ NOT MATCHED")

    fp = len(pred) - tp
    fn = len(truth) - tp

    prec = tp / (tp + fp) if (tp + fp) else 0
    rec = tp / (tp + fn) if (tp + fn) else 0
    f1 = 2 * prec * rec / (prec + rec) if (prec + rec) else 0

    return tp, fp, fn, prec, rec, f1, matched_details


###############################################################################
#                                      MAIN
###############################################################################

def main():
    ap = argparse.ArgumentParser(description="Evaluate long repeat predictions.")

    ap.add_argument("--truth", required=True)
    ap.add_argument("--pred", required=True)
    ap.add_argument("--overlap", type=float, default=0.8)

    ap.add_argument("--truth-format",
                    choices=["groundtruth", "minimap"],
                    default="groundtruth")

    ap.add_argument("--pred-format",
                    choices=["repeatfinder", "minimap", "mummer"],
                    default="repeatfinder")

    ap.add_argument("--minlen", type=int, default=15000)

    args = ap.parse_args()

    print("=" * 80)
    print("STEP 1: Loading TRUTH")
    print("=" * 80)

    if args.truth_format == "minimap":
        truth = read_minimap_repeats(args.truth, minlen=args.minlen)
    else:
        truth = read_truth_repeats(args.truth)

    print(f"\n[INFO] Truth count = {len(truth)}")

    print("\n" + "=" * 80)
    print("STEP 2: Loading PREDICTIONS")
    print("=" * 80)

    if args.pred_format == "minimap":
        pred = read_minimap_repeats(args.pred, minlen=args.minlen)
    elif args.pred_format == "mummer":
        pred = read_mummer_repeats(args.pred, minlen=args.minlen)
    else:
        pred = read_pred_repeats(args.pred)

    print(f"\n[INFO] Pred count = {len(pred)}")

    if not truth:
        print("[ERROR] No truth intervals loaded!")
        return
    if not pred:
        print("[ERROR] No predictions loaded!")
        return

    print("\n" + "=" * 80)
    print("STEP 3: Evaluating")
    print("=" * 80)

    tp, fp, fn, prec, rec, f1, matches = evaluate_truth_vs_pred(
        truth, pred, args.overlap
    )

    print("\n" + "=" * 80)
    print("FINAL RESULTS")
    print("=" * 80)
    print(f"TP = {tp}")
    print(f"FP = {fp}")
    print(f"FN = {fn}")
    print(f"Precision = {prec:.3f}")
    print(f"Recall    = {rec:.3f}")
    print(f"F1        = {f1:.3f}")

    if matches:
        print("\nMatched pairs:")
        for (ti, pi, t, p, s) in matches:
            print(f"  Truth #{ti+1} <-> Pred #{pi}, score={s:.3f}")


if __name__ == "__main__":
    main()
