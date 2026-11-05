#!/usr/bin/env python3
"""NOTE: the built-in `auto_min_abund` below returns 2 for almost every sample; run with --cutoffs (gmg-core's rule).

Read-set sketches (.gsk, abundance) -> QC, decontamination, distances, groups, ordination, NJ tree, label checks.

usage: 02_popdist.py SKETCH_DIR MANIFEST.tsv OUT_DIR [options]

Method (alignment-free, one pass over the sketches):
  1. per sample: drop hashes below an auto min-abundance (valley below the main coverage peak), which removes most
     sequencing-error k-mers;
  2. optional decontamination of flagged samples (--flag): keep a hash only if it is in a reference sketch of the
     target genus (--ref-hashes) or seen in >= --min-share unflagged samples (the cohort votes). Unflagged samples
     are left untouched, so clean samples are not subject to any ascertainment filter;
  3. containment C(a,b) = |A&B| / min(|A|,|B|); ANI ~ C**(1/k); mash distance from Jaccard also written;
  4. UPGMA groups at a distance cut, PCoA, DAPC-style (PCA -> LDA) on the hash presence matrix, NJ tree,
     per-sample label check (median ANI to each label) and nearest neighbours.
"""
import argparse, csv, os, struct, sys
import numpy as np
from scipy import sparse


def read_gsk(path):
    with open(path, "rb") as f:
        assert f.read(8) == b"GMGSKT01"
        k, scaled, seed, nl = struct.unpack("<IIQI", f.read(20))
        name = f.read(nl).decode()
        (n,) = struct.unpack("<Q", f.read(8))
        h = np.frombuffer(f.read(8 * n), dtype="<u8")
        has = f.read(1)[0]
        a = np.frombuffer(f.read(4 * n), dtype="<u4") if has else None
    return k, scaled, name, h, a


def auto_min_abund(a):
    """Valley between the error peak (abundance 1) and the main coverage peak (searched over the whole range)."""
    hist = np.bincount(np.minimum(a, 2000), minlength=2001)[1:].astype(float)  # index i = abundance i+1
    sm = np.convolve(hist, np.ones(5) / 5, mode="same")
    peak = 2 + int(np.argmax(sm[2:]))
    if sm[peak] < 3:
        return 2
    hi = max(2, peak // 2)
    return 1 + int(np.argmin(sm[1:hi + 1])) + 1


def nj(D, names):
    n = len(names)
    D = D.astype(float).copy()
    nodes = list(names)
    active = list(range(n))
    while len(active) > 2:
        m = len(active)
        sub = D[np.ix_(active, active)]
        r = sub.sum(1)
        Q = (m - 2) * sub - r[:, None] - r[None, :]
        np.fill_diagonal(Q, np.inf)
        i, j = np.unravel_index(np.argmin(Q), Q.shape)
        dij = sub[i, j]
        li = 0.5 * dij + (r[i] - r[j]) / (2 * (m - 2))
        lj = dij - li
        ai, aj = active[i], active[j]
        nodes[ai] = f"({nodes[ai]}:{max(li, 0):.6f},{nodes[aj]}:{max(lj, 0):.6f})"
        for t in active:
            if t not in (ai, aj):
                D[ai, t] = D[t, ai] = 0.5 * (D[ai, t] + D[aj, t] - dij)
        active.remove(aj)
    a, b = active
    return f"({nodes[a]}:{D[a, b] / 2:.6f},{nodes[b]}:{D[a, b] / 2:.6f});"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("sketch_dir"); ap.add_argument("manifest"); ap.add_argument("out")
    ap.add_argument("--exclude", default="", help="comma list of sample ids to drop entirely")
    ap.add_argument("--only-label", default="", help="keep only samples whose (relabelled) organism equals this")
    ap.add_argument("--relabel", default="", help="comma list id=Organism overriding the manifest label (e.g. confirmed swaps)")
    ap.add_argument("--cutoffs", default="", help="TSV with columns id and rust_min_abund (cutoff_probe output): use gmg-core's "
                    "auto_min_abund per sample instead of the built-in rule, which in practice returns 2 for almost every sample")
    ap.add_argument("--flag", default="", help="file or comma list of sample ids to decontaminate")
    ap.add_argument("--ref-hashes", default="", help=".npy of target-genus reference hashes (same k/scaled)")
    ap.add_argument("--min-share", type=int, default=3, help="decontam: keep hashes seen in this many unflagged samples")
    ap.add_argument("--min-hashes", type=int, default=2000, help="drop samples with fewer filtered hashes")
    ap.add_argument("--min-cov", type=float, default=4.0, help="drop samples with median k-mer coverage below this")
    ap.add_argument("--cut", type=float, default=0.02, help="UPGMA cut on (1-ANI)")
    ap.add_argument("--max-cols", type=int, default=200000, help="max hash columns used for PCA/DAPC")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    excl = {x for x in args.exclude.split(",") if x}
    flag = set()
    if args.flag:
        flag = set(open(args.flag).read().split()) if os.path.exists(args.flag) else set(args.flag.split(","))
    meta = {r["id"]: r for r in csv.DictReader(open(args.manifest), delimiter="\t")}
    cut = {}
    if args.cutoffs:
        cut = {r["id"]: int(r["rust_min_abund"]) for r in csv.DictReader(open(args.cutoffs), delimiter="\t")}
    for kv in filter(None, args.relabel.split(",")):
        i, lab = kv.split("=", 1)
        meta[i]["organism"] = lab

    ids, sets, qc, k, scaled = [], [], [], None, None
    for fn in sorted(os.listdir(args.sketch_dir)):
        if not fn.endswith(".gsk"):
            continue
        sid = fn[:-4]
        k, scaled, _, h, a = read_gsk(os.path.join(args.sketch_dir, fn))
        thr = cut.get(sid, auto_min_abund(a))
        keep = a >= thr
        med = float(np.median(a[keep])) if keep.any() else 0.0
        row = dict(id=sid, organism=meta.get(sid, {}).get("organism", ""), raw_hashes=len(h), min_abund=thr,
                   kept_hashes=int(keep.sum()), est_genome_mb=keep.sum() * scaled / 1e6, median_kmer_cov=med,
                   flagged=int(sid in flag), decontam_hashes=0, decontam_mb=0.0, status="ok")
        if sid in excl or (args.only_label and row["organism"] != args.only_label):
            row["status"] = "excluded"
        elif med < args.min_cov:
            row["status"] = "low_coverage"
        qc.append(row)
        if row["status"] == "ok":
            ids.append(sid); sets.append(h[keep])

    # decontamination of flagged samples
    if flag and args.ref_hashes:
        ref = np.load(args.ref_hashes)
        clean = [s for sid, s in zip(ids, sets) if sid not in flag]
        uh, uc = np.unique(np.concatenate(clean), return_counts=True)
        vote = uh[uc >= args.min_share]
        allowed = np.union1d(ref, vote)
        qmap = {r["id"]: r for r in qc}
        for i, sid in enumerate(ids):
            if sid in flag:
                new = sets[i][np.isin(sets[i], allowed)]
                qmap[sid]["decontam_hashes"] = len(new)
                qmap[sid]["decontam_mb"] = len(new) * scaled / 1e6
                sets[i] = new
    for r in qc:
        if r["status"] == "ok" and len(sets[ids.index(r["id"])]) < args.min_hashes:
            r["status"] = "too_few_hashes"
    keep_i = [i for i, sid in enumerate(ids) if next(r for r in qc if r["id"] == sid)["status"] == "ok"]
    ids = [ids[i] for i in keep_i]; sets = [sets[i] for i in keep_i]
    with open(f"{args.out}/qc.tsv", "w") as f:
        w = csv.DictWriter(f, fieldnames=list(qc[0]), delimiter="\t"); w.writeheader(); w.writerows(qc)
    n = len(ids)
    print(f"{n} samples kept of {len(qc)}", file=sys.stderr)
    labels = np.array([meta.get(s, {}).get("organism", "") for s in ids])

    # sparse presence matrix
    allh, cnt = np.unique(np.concatenate(sets), return_counts=True)
    sizes = np.array([len(s) for s in sets])
    rows = np.repeat(np.arange(n), sizes)
    cols = np.searchsorted(allh, np.concatenate(sets))
    S = sparse.csr_matrix((np.ones(len(rows), dtype=np.float32), (rows, cols)), shape=(n, len(allh)))
    inter = (S @ S.T).toarray().astype(float)
    sz = sizes.astype(float)
    C = inter / np.minimum(sz[:, None], sz[None, :])
    J = inter / (sz[:, None] + sz[None, :] - inter)
    ani = np.clip(C, 1e-9, 1) ** (1.0 / k)
    np.fill_diagonal(ani, 1.0)
    mash = np.where(J > 0, -np.log(2 * J / (1 + J)) / k, 1.0)
    np.fill_diagonal(mash, 0)
    for name, mat in (("containment", C), ("jaccard", J), ("ani_est", ani), ("mash_dist", mash)):
        np.savetxt(f"{args.out}/{name}.tsv", mat, delimiter="\t", fmt="%.5f", header="\t".join(ids), comments="")

    iu = np.triu_indices(n, 1)
    h, edges = np.histogram(ani[iu], bins=np.arange(0.90, 1.0001, 0.0025))
    with open(f"{args.out}/ani_hist.tsv", "w") as f:
        f.write("ani_lo\tpairs\n"); [f.write(f"{e:.4f}\t{c}\n") for e, c in zip(edges[:-1], h)]

    # groups
    from scipy.cluster.hierarchy import linkage, fcluster
    from scipy.spatial.distance import squareform
    Dm = 1 - ani; np.fill_diagonal(Dm, 0); Dm = (Dm + Dm.T) / 2
    grp = fcluster(linkage(squareform(Dm, checks=False), "average"), args.cut, "distance")
    with open(f"{args.out}/groups.tsv", "w") as f:
        f.write("id\torganism_label\tgroup\n")
        for s, l, g in zip(ids, labels, grp):
            f.write(f"{s}\t{l}\t{g}\n")

    # label check: median ANI to each label (self excluded) + nearest neighbours
    uniq = [l for l in np.unique(labels) if (labels == l).sum() >= 3]
    with open(f"{args.out}/label_check.tsv", "w") as f, open(f"{args.out}/neighbors.tsv", "w") as g2:
        f.write("id\tlabel\t" + "\t".join(f"median_ani_to_{l.replace(' ', '_')}" for l in uniq) + "\tclosest_label\tconflict\n")
        g2.write("id\tlabel\tnn1\tnn1_ani\tnn1_label\tnn2\tnn2_ani\tnn3\tnn3_ani\tnn4\tnn4_ani\tnn5\tnn5_ani\n")
        for i, s in enumerate(ids):
            med = []
            for l in uniq:
                m = (labels == l); m[i] = False
                med.append(float(np.median(ani[i, m])))
            best = uniq[int(np.argmax(med))]
            f.write(f"{s}\t{labels[i]}\t" + "\t".join(f"{x:.4f}" for x in med) + f"\t{best}\t{int(best != labels[i])}\n")
            o = [j for j in np.argsort(-ani[i]) if j != i][:5]
            g2.write(f"{s}\t{labels[i]}\t{ids[o[0]]}\t{ani[i, o[0]]:.4f}\t{labels[o[0]]}\t" +
                     "\t".join(f"{ids[j]}\t{ani[i, j]:.4f}" for j in o[1:]) + "\n")

    # PCoA
    Jn = np.eye(n) - 1 / n
    B = -0.5 * Jn @ (Dm ** 2) @ Jn
    w, v = np.linalg.eigh(B); o = np.argsort(w)[::-1][:10]
    pc = v[:, o] * np.sqrt(np.maximum(w[o], 0))
    with open(f"{args.out}/pcoa.tsv", "w") as f:
        f.write("id\t" + "\t".join(f"PC{i+1}" for i in range(10)) + "\n")
        for s, row in zip(ids, pc):
            f.write(s + "\t" + "\t".join(f"{x:.5f}" for x in row) + "\n")
    with open(f"{args.out}/pcoa_var.txt", "w") as f:
        f.write(" ".join(f"{x:.4f}" for x in np.maximum(w[o], 0) / np.maximum(w, 0).sum()))

    # DAPC-style: PCA on presence/absence of shared hashes, LDA on UPGMA groups with >= 3 members
    try:
        from sklearn.decomposition import PCA
        from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
        cand = np.where((cnt >= 2) & (cnt <= n - 1))[0]
        if len(cand) > args.max_cols:
            cand = np.sort(np.random.default_rng(0).choice(cand, args.max_cols, replace=False))
        X = S[:, cand].toarray(); X -= X.mean(0)
        P = PCA(n_components=min(30, n - 1), svd_solver="randomized", random_state=0).fit_transform(X)
        del X
        vals, counts = np.unique(grp, return_counts=True)
        big = np.isin(grp, vals[counts >= 3])
        if len(np.unique(grp[big])) >= 2:
            ld = LinearDiscriminantAnalysis().fit(P[big], grp[big]).transform(P)
            with open(f"{args.out}/dapc.tsv", "w") as f:
                f.write("id\tgroup\tlabel\t" + "\t".join(f"LD{i+1}" for i in range(ld.shape[1])) + "\n")
                for s, g, l, row in zip(ids, grp, labels, ld):
                    f.write(f"{s}\t{g}\t{l}\t" + "\t".join(f"{x:.4f}" for x in row) + "\n")
        with open(f"{args.out}/pca_hash.tsv", "w") as f:
            f.write("id\t" + "\t".join(f"PC{i+1}" for i in range(5)) + "\n")
            for s, row in zip(ids, P):
                f.write(s + "\t" + "\t".join(f"{x:.3f}" for x in row[:5]) + "\n")
        print(f"PCA/DAPC on {len(cand)} hash columns", file=sys.stderr)
    except Exception as e:
        print("DAPC skipped:", e, file=sys.stderr)

    nwk = nj(mash, ids)
    open(f"{args.out}/nj.nwk", "w").write(nwk + "\n")
    try:  # midpoint-rooted species monophyly check + figure
        from io import StringIO
        from Bio import Phylo
        import matplotlib; matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        t = Phylo.read(StringIO(nwk), "newick"); t.root_at_midpoint()
        lab = dict(zip(ids, labels)); res = []
        for l in np.unique(labels):
            names = [s for s in ids if lab[s] == l]
            if len(names) < 2:
                continue
            tips = [c for c in t.get_terminals() if c.name in set(names)]
            mono = t.is_monophyletic(tips)
            res.append(f"{l}\tn={len(names)}\tmonophyletic={mono is not False}")
        open(f"{args.out}/tree_check.txt", "w").write("\n".join(res) + "\n")
        pal = {l: c for l, c in zip(np.unique(labels), ["#d95f02", "#1b9e77", "#7570b3", "#e7298a"])}
        for c in t.get_terminals():
            c.color = pal.get(lab[c.name], "gray")
        fig, ax = plt.subplots(figsize=(8, 12))
        Phylo.draw(t, axes=ax, label_func=lambda c: None, do_show=False)
        ax.set_title("NJ on mash distance (midpoint root); orange=immitis green=posadasii")
        fig.savefig(f"{args.out}/nj_tree.png", dpi=90)
    except Exception as e:
        print("tree check skipped:", e, file=sys.stderr)
    print("wrote", args.out, file=sys.stderr)


if __name__ == "__main__":
    main()
