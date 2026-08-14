#!/usr/bin/env python3
"""
Genomovar analysis within the major Sulfitobacter clades.

Addresses reviewer request: test whether each major clade subdivides into
genomovars across the reported 99.2-99.8% within-species ANI gap
(Rodriguez-R et al.), and whether such genomovars are fixed for the
ILS-diagnostic (ancestral) sequence variants.

Inputs (all in this directory except the alignment; see README.md)
------
ANI      : Sulfito_ANI.txt                     (fastANI, 74 genomes)
Strains  : non_redundant_gnm.txt                (45 non-duplicate, clade)
Alignment: focal_clades_concate_coreLCB.fasta
           (34 non-duplicate members of C1/C2/C9/C12, 2,803,078 bp core LCB;
           not included in this repository, see README.md for how to obtain it)

Analyses are restricted to the 45 non-duplicate strains throughout.

Usage: python3 genomovar_analysis.py
"""
import csv
import itertools
import json
import os
import sys
from collections import Counter, defaultdict

import numpy as np
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform

HERE = os.path.dirname(os.path.abspath(__file__))
ANI_FILE = os.path.join(HERE, "Sulfito_ANI.txt")
STRAIN_FILE = os.path.join(HERE, "non_redundant_gnm.txt")
ALN_FILE = os.path.join(HERE, "focal_clades_concate_coreLCB.fasta")
META_FILE = os.path.join(HERE, "metadata.tsv")
OUTDIR = HERE

# Genomovar boundary: Rodriguez-R et al. report a depletion of within-species
# ANI values between 99.2% and 99.8%; 99.5% is the operational midpoint.
GAP_LO, GAP_HI, GV_THRESHOLD = 99.2, 99.8, 99.5
LINKAGES = ("average", "single", "complete")


# --------------------------------------------------------------------------- #
# I/O
# --------------------------------------------------------------------------- #
def load_strains():
    """Return {strain: clade} for the 45 non-duplicate genomes."""
    out = {}
    with open(STRAIN_FILE) as fh:
        rdr = csv.reader(fh, delimiter="\t")
        next(rdr)  # header: Clade / strain
        for row in rdr:
            if len(row) >= 2 and row[1].strip():
                out[row[1].strip()] = row[0].strip()
    return out


def load_metadata():
    """Return {strain: (clade, ecological source, geographical source)}."""
    out = {}
    with open(META_FILE) as fh:
        rdr = csv.reader(fh, delimiter="\t")
        next(rdr)
        for row in rdr:
            if len(row) >= 5:
                out[row[1].strip()] = (row[0].strip(), row[3].strip(), row[4].strip())
    return out


def load_ani():
    """Return symmetrised ANI lookup. fastANI is directional; we average both
    directions so that the matrix used for clustering is symmetric."""
    raw = {}
    with open(ANI_FILE) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 3:
                continue
            a = p[0].split("/")[-1].replace(".fasta", "")
            b = p[1].split("/")[-1].replace(".fasta", "")
            raw[(a, b)] = float(p[2])

    def get(a, b):
        vals = [v for v in (raw.get((a, b)), raw.get((b, a))) if v is not None]
        return float(np.mean(vals)) if vals else None

    return raw, get


def load_alignment(keep):
    """Read the core-LCB concatenation, returning {strain: np.uint8 array}."""
    seqs, name, buf = {}, None, []
    with open(ALN_FILE) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None and name in keep:
                    seqs[name] = np.frombuffer(
                        "".join(buf).upper().encode(), dtype=np.uint8
                    )
                name, buf = line[1:].strip().split()[0], []
            else:
                buf.append(line.strip())
    if name is not None and name in keep:
        seqs[name] = np.frombuffer("".join(buf).upper().encode(), dtype=np.uint8)
    return seqs


# --------------------------------------------------------------------------- #
# Clustering
# --------------------------------------------------------------------------- #
def ani_matrix(members, get):
    """Symmetric ANI matrix (percent) for an ordered list of strains."""
    n = len(members)
    mat = np.full((n, n), np.nan)
    for i in range(n):
        mat[i, i] = 100.0
    for i, j in itertools.combinations(range(n), 2):
        v = get(members[i], members[j])
        if v is not None:
            mat[i, j] = mat[j, i] = v
    return mat


def cluster_at(members, get, threshold, method="average"):
    """Cluster genomes so that members of a cluster are >= threshold ANI.
    Distance = 100 - ANI; cut height = 100 - threshold."""
    if len(members) == 1:
        return {members[0]: 1}
    mat = ani_matrix(members, get)
    dist = 100.0 - mat
    np.fill_diagonal(dist, 0.0)
    # fastANI omits pairs below ~80% identity; such pairs are maximally distant.
    dist[np.isnan(dist)] = 100.0
    dist = (dist + dist.T) / 2.0
    Z = linkage(squareform(dist, checks=False), method=method)
    labels = fcluster(Z, t=100.0 - threshold, criterion="distance")
    # Relabel so cluster 1 is the largest, for stable presentation.
    order = [c for c, _ in Counter(labels).most_common()]
    remap = {c: k + 1 for k, c in enumerate(order)}
    return {m: remap[l] for m, l in zip(members, labels)}


def within_pairs(members, get):
    """Sorted list of within-group pairwise ANI values."""
    vals = [get(a, b) for a, b in itertools.combinations(sorted(members), 2)]
    return sorted(v for v in vals if v is not None)


# --------------------------------------------------------------------------- #
# Sequence-level tests
# --------------------------------------------------------------------------- #
BASES = np.frombuffer(b"ACGT", dtype=np.uint8)


def unambiguous_mask(seqs, members):
    """Boolean mask of alignment columns where every member has an A/C/G/T."""
    arr = np.vstack([seqs[m] for m in members])
    is_base = np.isin(arr, BASES)          # (n_rows, n_cols)
    mask = is_base.all(axis=0)             # column usable iff all rows are ACGT
    return arr, mask


def nucleotide_diversity(arr, mask):
    """Mean pairwise difference per usable site (pi) for the rows of arr."""
    sub = arr[:, mask]
    n = sub.shape[0]
    if n < 2 or sub.shape[1] == 0:
        return 0.0, 0
    tot = 0
    for i, j in itertools.combinations(range(n), 2):
        tot += int(np.count_nonzero(sub[i] != sub[j]))
    npairs = n * (n - 1) // 2
    return tot / npairs / sub.shape[1], sub.shape[1]


def polymorphic_columns(arr, mask):
    """Indices (into the full alignment) of columns variable among rows of arr."""
    sub = arr[:, mask]
    idx = np.flatnonzero(mask)
    var = (sub != sub[0]).any(axis=0)
    return idx[var], sub[:, var]


def genomovar_partition_test(members, gv, sub_cols, col_idx):
    """Reviewer's prediction: if genomovars are the unit that is fixed for
    ancestral (ILS) variants, then within-clade polymorphism should be
    *between* genomovars (each genomovar internally fixed), not within them.

    Returns counts of within-clade polymorphic sites that are
      - 'between'  : every genomovar internally fixed (differential fixation)
      - 'within'   : at least one genomovar still polymorphic
    """
    groups = defaultdict(list)
    for k, m in enumerate(members):
        groups[gv[m]].append(k)
    multi = [rows for rows in groups.values() if len(rows) > 1]
    if len(groups) < 2:
        return dict(between=0, within=int(sub_cols.shape[1]), n_gv=len(groups))
    between = 0
    within = 0
    for c in range(sub_cols.shape[1]):
        col = sub_cols[:, c]
        internally_fixed = all(len(set(col[rows].tolist())) == 1 for rows in multi)
        if internally_fixed:
            between += 1
        else:
            within += 1
    return dict(between=between, within=within, n_gv=len(groups))


def informative_partition_test(members, gv, sub_cols, n_perm=1000, seed=0):
    """The raw 'between genomovar' count is inflated at high ANI thresholds
    because singleton genomovars are fixed by definition. Restrict the test to
    strains sitting in multi-member genomovars (where 'internally fixed' is a
    real constraint), and compare the observed count of differentially fixed
    sites to a permutation null that preserves genomovar sizes.
    """
    groups = defaultdict(list)
    for k, m in enumerate(members):
        groups[gv[m]].append(k)
    multi = {g: rows for g, rows in groups.items() if len(rows) > 1}
    rows_used = sorted(r for rows in multi.values() for r in rows)
    if len(multi) < 2 or len(rows_used) < 4:
        return dict(testable=False, n_gv_multi=len(multi), n_strains=len(rows_used),
                    observed=None, null_mean=None, null_sd=None, p=None, n_sites=0)
    sub = sub_cols[rows_used]
    remap = {r: i for i, r in enumerate(rows_used)}
    blocks = [[remap[r] for r in rows] for rows in multi.values()]
    var = (sub != sub[0]).any(axis=0)
    sub = sub[:, var]

    def count_fixed(assign_blocks):
        fixed = np.ones(sub.shape[1], dtype=bool)
        for rows in assign_blocks:
            first = sub[rows[0]]
            same = np.ones(sub.shape[1], dtype=bool)
            for r in rows[1:]:
                same &= sub[r] == first
            fixed &= same
        return int(np.count_nonzero(fixed))

    obs = count_fixed(blocks)
    rng = np.random.default_rng(seed)
    sizes = [len(b) for b in blocks]
    null = np.empty(n_perm, dtype=int)
    idx = np.arange(sub.shape[0])
    for p in range(n_perm):
        perm = rng.permutation(idx)
        cut, rand_blocks = 0, []
        for s in sizes:
            rand_blocks.append(list(perm[cut:cut + s]))
            cut += s
        null[p] = count_fixed(rand_blocks)
    pval = (1 + int(np.count_nonzero(null >= obs))) / (n_perm + 1)
    return dict(testable=True, n_gv_multi=len(multi), n_strains=len(rows_used),
                observed=obs, null_mean=float(null.mean()),
                null_sd=float(null.std()), p=pval, n_sites=int(sub.shape[1]))


def trans_clade_shared(seqs, clade_members, focal, gv):
    """Sites polymorphic within the focal clade whose minor allele is also
    present in another clade -> candidate retained ancestral polymorphism.
    Reports how those sites partition with respect to genomovars."""
    members = clade_members[focal]
    others = [m for c, ms in clade_members.items() if c != focal for m in ms]
    all_m = members + others
    arr_all, mask_all = unambiguous_mask(seqs, all_m)
    nf = len(members)
    sub = arr_all[:, mask_all]
    idx = np.flatnonzero(mask_all)
    foc = sub[:nf]
    oth = sub[nf:]
    var = (foc != foc[0]).any(axis=0)
    foc_v, oth_v, idx_v = foc[:, var], oth[:, var], idx[var]
    shared_cols, shared_idx = [], []
    for c in range(foc_v.shape[1]):
        alleles = set(foc_v[:, c].tolist())
        if len(alleles) != 2:
            continue
        counts = Counter(foc_v[:, c].tolist())
        minor = counts.most_common()[-1][0]
        if np.any(oth_v[:, c] == minor):
            shared_cols.append(c)
            shared_idx.append(int(idx_v[c]))
    if not shared_cols:
        return dict(n_shared=0, between=0, within=0, positions=[])
    part = genomovar_partition_test(
        members, gv, foc_v[:, shared_cols], np.array(shared_idx)
    )
    return dict(
        n_shared=len(shared_cols),
        between=part["between"],
        within=part["within"],
        positions=shared_idx,
    )


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #
def main():
    strains = load_strains()
    meta = load_metadata()
    raw, get = load_ani()
    print("Non-duplicate strains: %d" % len(strains))
    assert len(strains) == 45, len(strains)
    for s in strains:
        assert meta[s][0] == strains[s], (s, meta[s][0], strains[s])
    print("Clade labels consistent with metadata.tsv: yes")

    clades = defaultdict(list)
    for s, c in strains.items():
        clades[c].append(s)
    for c in clades:
        clades[c].sort()
    multi = sorted(
        [c for c in clades if len(clades[c]) > 1], key=lambda x: -len(clades[x])
    )
    print("Multi-member clades: %s" % ", ".join("%s(n=%d)" % (c, len(clades[c])) for c in multi))

    # ---- 1. within-clade ANI distribution -------------------------------- #
    rows = []
    pooled = []
    for c in multi:
        v = within_pairs(clades[c], get)
        pooled.extend(v)
        rows.append(
            dict(
                clade=c,
                n_genomes=len(clades[c]),
                n_pairs=len(v),
                min=min(v),
                max=max(v),
                mean=float(np.mean(v)),
                median=float(np.median(v)),
                n_in_gap=sum(1 for x in v if GAP_LO <= x <= GAP_HI),
            )
        )
    print("\n--- within-clade pairwise ANI (45 non-duplicate strains) ---")
    print("%-5s %5s %6s %8s %8s %8s %8s %10s" % ("clade","n","pairs","min","max","mean","median","in 99.2-99.8"))
    for r in rows:
        print("%-5s %5d %6d %8.2f %8.2f %8.2f %8.2f %6d (%.0f%%)" % (
            r["clade"], r["n_genomes"], r["n_pairs"], r["min"], r["max"],
            r["mean"], r["median"], r["n_in_gap"], 100*r["n_in_gap"]/r["n_pairs"]))
    pooled = np.array(pooled)
    n_gap = int(np.sum((pooled >= GAP_LO) & (pooled <= GAP_HI)))
    print("POOLED n=%d  min=%.2f max=%.2f  in gap=%d (%.1f%%)" % (
        len(pooled), pooled.min(), pooled.max(), n_gap, 100*n_gap/len(pooled)))

    with open(os.path.join(OUTDIR, "table_within_clade_ANI.tsv"), "w") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    np.savetxt(os.path.join(OUTDIR, "pooled_within_clade_ANI.txt"), pooled, fmt="%.6f")

    # ---- 2. genomovar delineation, threshold/linkage sensitivity --------- #
    sens = []
    for method in LINKAGES:
        for thr in (99.2, 99.5, 99.8):
            for c in multi:
                gvv = cluster_at(clades[c], get, thr, method)
                sens.append(dict(linkage=method, threshold=thr, clade=c,
                                 n_genomes=len(clades[c]),
                                 n_genomovars=len(set(gvv.values())),
                                 sizes=",".join(str(n) for _, n in Counter(gvv.values()).most_common())))
    with open(os.path.join(OUTDIR, "table_genomovar_sensitivity.tsv"), "w") as fh:
        w = csv.DictWriter(fh, fieldnames=list(sens[0]), delimiter="\t")
        w.writeheader()
        for r in sens:
            w.writerow(r)
    print("\n--- genomovar counts by threshold (average linkage) ---")
    print("%-5s %5s %8s %8s %8s" % ("clade","n","99.2%","99.5%","99.8%"))
    for c in multi:
        cells = []
        for thr in (99.2, 99.5, 99.8):
            k = [r for r in sens if r["linkage"]=="average" and r["threshold"]==thr and r["clade"]==c][0]
            cells.append("%d (%s)" % (k["n_genomovars"], k["sizes"]))
        print("%-5s %5d %8s %8s %8s" % (c, len(clades[c]), cells[0], cells[1], cells[2]))

    # ---- 3. genomovar assignment at the operational threshold ------------ #
    gv_all = {}
    for c in multi:
        gvv = cluster_at(clades[c], get, GV_THRESHOLD, "average")
        for m, g in gvv.items():
            gv_all[m] = "%s.gv%d" % (c, g)
    for c in clades:
        if len(clades[c]) == 1:
            gv_all[clades[c][0]] = "%s.gv1" % c
    with open(os.path.join(OUTDIR, "table_genomovar_assignment.tsv"), "w") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["clade","strain","genomovar_99.5","ecological_source","geographical_source"])
        for c in sorted(clades, key=lambda x: (-len(clades[x]), x)):
            for s in clades[c]:
                w.writerow([c, s, gv_all[s], meta[s][1], meta[s][2]])
    print("\nWrote genomovar assignments for %d strains" % len(gv_all))

    # ---- 4. sequence-level test on the four focal clades ----------------- #
    focal = ["C1", "C2", "C9", "C12"]
    keep = set(s for c in focal for s in clades[c])
    seqs = load_alignment(keep)
    print("\n--- core-LCB alignment ---")
    print("strains loaded: %d ; expected: %d" % (len(seqs), len(keep)))
    assert set(seqs) == keep, sorted(keep - set(seqs))
    lens = set(int(v.size) for v in seqs.values())
    assert len(lens) == 1, lens
    print("alignment length: %d bp (equal across strains)" % lens.pop())

    clade_members = {c: clades[c] for c in focal}
    seq_rows = []
    for c in focal:
        members = clade_members[c]
        arr, mask = unambiguous_mask(seqs, members)
        pi, nsites = nucleotide_diversity(arr, mask)
        pos, sub = polymorphic_columns(arr, mask)
        gvv = cluster_at(members, get, GV_THRESHOLD, "average")
        part = genomovar_partition_test(members, gvv, sub, pos)
        tcs = trans_clade_shared(seqs, clade_members, c, gvv)
        seq_rows.append(dict(
            clade=c, n_genomes=len(members), usable_sites=nsites, pi=pi,
            n_polymorphic=int(sub.shape[1]), n_genomovars=part["n_gv"],
            poly_between_gv=part["between"], poly_within_gv=part["within"],
            trans_clade_shared=tcs["n_shared"],
            tcs_between_gv=tcs["between"], tcs_within_gv=tcs["within"]))

    print("\n--- within-clade polymorphism vs genomovar structure (99.5%, average) ---")
    hdr = ("clade","n","pi","poly sites","n_gv","between gv","within gv","transclade","tc betw","tc with")
    print("%-5s %3s %10s %11s %5s %11s %10s %11s %8s %8s" % hdr)
    for r in seq_rows:
        print("%-5s %3d %10.2e %11d %5d %11d %10d %11d %8d %8d" % (
            r["clade"], r["n_genomes"], r["pi"], r["n_polymorphic"], r["n_genomovars"],
            r["poly_between_gv"], r["poly_within_gv"], r["trans_clade_shared"],
            r["tcs_between_gv"], r["tcs_within_gv"]))

    with open(os.path.join(OUTDIR, "table_sequence_level.tsv"), "w") as fh:
        w = csv.DictWriter(fh, fieldnames=list(seq_rows[0]), delimiter="\t")
        w.writeheader()
        for r in seq_rows:
            w.writerow(r)
    print("\nWrote table_sequence_level.tsv")

    # ---- 5. threshold sweep for the partition test ----------------------- #
    print("\n--- partition test across thresholds (C1 and C12, average linkage) ---")
    print("%-5s %6s %5s %11s %11s %9s %s" % (
        "clade","thresh","n_gv","between gv","within gv","% between","note"))
    sweep = []
    for c in ("C1", "C12", "C2"):
        members = clade_members[c]
        arr, mask = unambiguous_mask(seqs, members)
        pos, sub = polymorphic_columns(arr, mask)
        for thr in (99.2, 99.5, 99.8, 99.9):
            gvv = cluster_at(members, get, thr, "average")
            part = genomovar_partition_test(members, gvv, sub, pos)
            tot = part["between"] + part["within"]
            pct = 100.0 * part["between"] / tot if tot else 0.0
            sizes = ",".join(str(n) for _, n in Counter(gvv.values()).most_common())
            print("%-5s %6.1f %5d %11d %11d %8.1f%% %s" % (
                c, thr, part["n_gv"], part["between"], part["within"], pct, sizes))
            sweep.append(dict(clade=c, threshold=thr, n_genomovars=part["n_gv"],
                              poly_between_gv=part["between"],
                              poly_within_gv=part["within"],
                              pct_between=round(pct, 2), genomovar_sizes=sizes))
    with open(os.path.join(OUTDIR, "table_partition_sweep.tsv"), "w") as fh:
        w = csv.DictWriter(fh, fieldnames=list(sweep[0]), delimiter="\t")
        w.writeheader()
        for r in sweep:
            w.writerow(r)

    # ---- 4b. how completely is ancestral variation sorted within clades? -- #
    pooled_members = [s for c in focal for s in clades[c]]
    arr_all, mask_all = unambiguous_mask(seqs, pooled_members)
    sub_all = arr_all[:, mask_all]
    n_alleles = np.zeros(sub_all.shape[1], dtype=np.int8)
    for b in BASES:
        n_alleles += (sub_all == b).any(axis=0)
    bi = n_alleles == 2
    sub_bi = sub_all[:, bi]
    row_of = {}
    start = 0
    for c in focal:
        row_of[c] = list(range(start, start + len(clades[c])))
        start += len(clades[c])
    fixed = np.ones(sub_bi.shape[1], dtype=bool)
    for c in focal:
        rows = row_of[c]
        same = np.ones(sub_bi.shape[1], dtype=bool)
        for r in rows[1:]:
            same &= sub_bi[r] == sub_bi[rows[0]]
        fixed &= same
    n_fixed = int(fixed.sum())
    print("\n--- completeness of sorting within clades ---")
    print("usable sites (all %d strains ACGT): %d" % (len(pooled_members), sub_all.shape[1]))
    print("biallelic sites: %d" % int(bi.sum()))
    print("fixed within all four focal clades: %d (%.1f%%)"
          % (n_fixed, 100.0 * n_fixed / int(bi.sum())))
    with open(os.path.join(OUTDIR, "table_sorting_completeness.tsv"), "w") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["n_strains","usable_sites","biallelic_sites",
                    "fixed_within_all_clades","percent_fixed"])
        w.writerow([len(pooled_members), sub_all.shape[1], int(bi.sum()), n_fixed,
                    round(100.0 * n_fixed / int(bi.sum()), 2)])

    # ---- 5b. permutation test on multi-member genomovars only ------------ #
    print("\n--- permutation test, multi-member genomovars only (1,000 permutations) ---")
    print("%-5s %6s %6s %8s %9s %11s %9s %8s" % (
        "clade","thresh","n_gv*","strains","sites","observed","null mean","P"))
    perm_rows = []
    for c in ("C1", "C12"):
        members = clade_members[c]
        arr, mask = unambiguous_mask(seqs, members)
        pos, sub = polymorphic_columns(arr, mask)
        for thr in (99.5, 99.8, 99.9):
            gvv = cluster_at(members, get, thr, "average")
            res = informative_partition_test(members, gvv, sub)
            if res["testable"]:
                print("%-5s %6.1f %6d %8d %9d %11d %9.1f %8.3f" % (
                    c, thr, res["n_gv_multi"], res["n_strains"], res["n_sites"],
                    res["observed"], res["null_mean"], res["p"]))
            else:
                print("%-5s %6.1f %6d %8d %9s %11s %9s %8s   not testable" % (
                    c, thr, res["n_gv_multi"], res["n_strains"], "-","-","-","-"))
            perm_rows.append(dict(clade=c, threshold=thr, **res))
    fields = ["clade","threshold","testable","n_gv_multi","n_strains","n_sites",
              "observed","null_mean","null_sd","p"]
    with open(os.path.join(OUTDIR, "table_permutation_test.tsv"), "w") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t",
                           extrasaction="ignore")
        w.writeheader()
        for r in perm_rows:
            w.writerow(r)

    # ---- 6. is C2's signal driven by one divergent strain? --------------- #
    print("\n--- C2 pairwise ANI (n=3) ---")
    for a, b in itertools.combinations(clades["C2"], 2):
        print("  %-10s %-10s %.2f" % (a, b, get(a, b)))
    print("\nWrote all tables to %s" % OUTDIR)
    return strains, meta, get, clades, multi, gv_all


if __name__ == "__main__":
    main()
