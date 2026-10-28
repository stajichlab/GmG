#!/usr/bin/env python3
"""Build docs/report-mockup/mockup.html from real `gmg run` output (docs/report-mockup/sample/).

This is a DESIGN MOCKUP (issue #22). The production engine (#21) will be written in Rust and must reproduce this
layout from results.jsonl/TSV. Everything shown is real data from a run of 148 inputs against the BFD index
(145 NCBI genomes not in BFD, two simulated read sets, one empty file). The only derived extras are the ANI matrix
and the abundance histogram, computed from sketches (the engine will emit both).

usage: make_mockup.py            (writes mockup.html next to this file)
"""
import csv, collections, gzip, html, json, math, os

HERE = os.path.dirname(os.path.abspath(__file__))
SAMPLE = os.path.join(HERE, "sample")
RANKS = ["phylum", "class", "order", "family", "genus", "species"]
TH = {"min_ani": [0.7109, 0.7109, 0.7109, 0.7715, 0.7765, 0.9938], "min_cons": [0.5, 0.5, 0.6, 0.9, 0.5, 0.9]}  # known
esc = html.escape

# Phylum colours follow the ENTITY (the phylum), never its rank in a run: the 8 phyla with the most BFD genomes
# take the 8 validated categorical slots in order; every other phylum shares one neutral grey (named in tooltips/tables).
PHYLUM_SLOTS = ["Ascomycota", "Basidiomycota", "Mucoromycota", "Kickxellomycota", "Microsporidia",
                "Mortierellomycota", "Chytridiomycota", "Glomeromycota"]


def ph_color(name):
    if name in PHYLUM_SLOTS:
        return f"var(--ph-{PHYLUM_SLOTS.index(name) + 1})"
    return "var(--nomatch)" if name in ("", "(not called)") else "var(--ph-other)"


def ph_order(name):
    return PHYLUM_SLOTS.index(name) if name in PHYLUM_SLOTS else (len(PHYLUM_SLOTS) if name not in ("", "(not called)") else len(PHYLUM_SLOTS) + 1)


def ph_dot(name):
    return f'<span class="dot" style="background:{ph_color(name)}" title="{esc(name or "phylum not called")}"></span>'


def read_tsv(name):
    with open(os.path.join(SAMPLE, name)) as f:
        return list(csv.DictReader(f, delimiter="\t"))


rows = read_tsv("results.tsv")
hits_rows = read_tsv("hits.tsv")
run = json.load(open(os.path.join(SAMPLE, "run.json")))
hist = json.load(open(os.path.join(SAMPLE, "abundance_hist.json")))
mat = json.load(open(os.path.join(SAMPLE, "ani_matrix.json")))
hits = collections.defaultdict(list)
for h in hits_rows:
    hits[h["query_id"]].append(h)
# cluster ids: the production index will carry them; here they come from the BFD cluster table (95% ANI)
cluster = {}
cl_path = os.path.join(HERE, "../../data/clusters/bfd_ani0.95.tsv.gz")
if os.path.exists(cl_path):
    with gzip.open(cl_path, "rt") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            cluster[r["id"]] = r["cluster"]
for r in rows:
    r["cluster"] = cluster.get(r["top_hit"], r.get("cluster", ""))
byid = {r["query_id"]: r for r in rows}


def f(x, d=4):
    try:
        return f"{float(x):.{d}f}"
    except (TypeError, ValueError):
        return ""


def short(s, n=26):
    return s if len(s) <= n else s[: n - 1] + "…"


# ------------------------------------------------------------------ small SVG helpers
def svg_text(x, y, s, cls="t-sec", anchor="start", extra=""):
    return f'<text x="{x:.1f}" y="{y:.1f}" class="{cls}" text-anchor="{anchor}" {extra}>{esc(str(s))}</text>'


def compact(n):
    n = float(n)
    if n >= 1e6:
        return f"{n/1e6:.1f}M".replace(".0M", "M")
    if n >= 1e3:
        return f"{n/1e3:.1f}K".replace(".0K", "K")
    return f"{n:g}"


def nice_ticks(vmax, n=4):
    if vmax <= 0:
        return [0, 1]
    raw = vmax / n
    mag = 10 ** math.floor(math.log10(raw))
    step = min((s * mag for s in (1, 2, 2.5, 5, 10) if s * mag >= raw), default=mag)
    return [i * step for i in range(int(vmax / step + 1.0001) + 1)]


# ------------------------------------------------------------------ run-level numbers
n_in = len(rows)
status_ct = collections.Counter(r["status"] for r in rows)
depth_ct = collections.Counter()
for r in rows:
    if r["status"] == "error":
        continue
    b = r["best_rank"]
    depth_ct[{"species": "species", "genus": "genus", "family": "family"}.get(b, "higher" if b else "none")] += 1
classified = depth_ct["species"] + depth_ct["genus"]
ok_total = n_in - status_ct["error"]

DEPTH = [("species", "Species", "var(--ord-1)", "#ffffff"), ("genus", "Genus", "var(--ord-2)", "#ffffff"),
         ("family", "Family", "var(--ord-3)", "#ffffff"), ("higher", "Order or higher", "var(--ord-4)", "#0b0b0b"),
         ("none", "No match", "var(--nomatch)", "var(--text-primary)")]


def depth_chart():
    W, H, x0, x1, y, h = 880, 78, 0, 880, 14, 30
    total = sum(depth_ct[k] for k, *_ in DEPTH)
    out, x, gap = [], x0, 2
    segs = []
    for k, label, fill, ink in DEPTH:
        n = depth_ct[k]
        if not n:
            continue
        w = (x1 - x0) * n / total
        tip = f"{label}: {n} of {total} inputs ({100*n/total:.0f}%)"
        out.append(f'<rect class="tip" tabindex="0" data-tip="{esc(tip)}" x="{x:.1f}" y="{y}" width="{max(w-gap,1):.1f}" height="{h}" rx="4" fill="{fill}"/>')
        if w > 70:  # label only when it fits with padding
            out.append(svg_text(x + w / 2 - gap / 2, y + h / 2 + 4.5, f"{n}", "t-fill", "middle", f'style="fill:{ink}"'))
        segs.append((k, label, n, fill))
        x += w
    # legend (direct labels would collide on the thin segments)
    lx = 0
    for k, label, n, fill in segs:
        out.append(f'<rect x="{lx}" y="{y+h+16}" width="10" height="10" rx="2" fill="{fill}"/>')
        txt = f"{label} {n}"
        out.append(svg_text(lx + 15, y + h + 25.5, txt, "t-sec"))
        lx += 15 + len(txt) * 6.6 + 18
    table = "<table><thead><tr><th>Deepest rank called</th><th>Inputs</th><th>Share</th></tr></thead><tbody>" + "".join(
        f"<tr><td>{esc(l)}</td><td>{n}</td><td>{100*n/total:.1f}%</td></tr>" for k, l, n, _ in segs) + "</tbody></table>"
    return f'<svg viewBox="0 0 {W} {H+6}" class="chart" role="img" aria-label="Deepest rank called per input">{"".join(out)}</svg>', table


def hbars(items, title_col, width=420, row_h=26, label_w=150, max_label=22, unit="", colors=None):
    """Single-series horizontal bars, value at the tip. items = [(label, value)]."""
    if not items:
        return "", ""
    vmax = max(v for _, v in items)
    plot_w = width - label_w - 44
    H = row_h * len(items) + 6
    out = []
    for i, (lab, v) in enumerate(items):
        y = i * row_h + 4
        bw = max(plot_w * v / vmax, 8)
        tip = f"{lab}: {v}{unit}"
        out.append(svg_text(label_w - 8, y + 14, short(lab, max_label), "t-sec", "end"))
        # 18px bar, 4px rounded data end, square at the baseline
        out.append(f'<path class="tip bar" tabindex="0" data-tip="{esc(tip)}" d="M{label_w},{y+3} h{bw-4:.1f} a4,4 0 0 1 4,4 v10 a4,4 0 0 1 -4,4 h-{bw-4:.1f} z" fill="{colors[i] if colors else 'var(--series-1)'}"/>')
        out.append(f'<rect class="tip hit" tabindex="-1" data-tip="{esc(tip)}" x="{label_w}" y="{y-1}" width="{plot_w+40}" height="{row_h-2}" fill="transparent"/>')
        out.append(svg_text(label_w + bw + 6, y + 14.5, f"{v}", "t-pri"))
    out.insert(0, f'<line x1="{label_w}" x2="{label_w}" y1="0" y2="{H}" class="axis"/>')
    table = f"<table><thead><tr><th>{esc(title_col)}</th><th>Inputs</th></tr></thead><tbody>" + "".join(
        f"<tr><td>{esc(l)}</td><td>{v}</td></tr>" for l, v in items) + "</tbody></table>"
    return f'<svg viewBox="0 0 {width} {H}" class="chart" role="img" aria-label="{esc(title_col)} counts">{"".join(out)}</svg>', table


def genus_by_phylum(items, width=420, row_h=24, head_h=26, label_w=130):
    """items = [(genus, phylum, n)]. Bars grouped under a phylum heading, coloured by phylum (heading + legend + tooltip + table
    all name the phylum, so identity never rests on colour alone)."""
    groups = collections.OrderedDict()
    for g, ph, n in sorted(items, key=lambda t: (ph_order(t[1]), -t[2], t[0])):
        groups.setdefault(ph or "(not called)", []).append((g, n))
    vmax = max(n for _, _, n in items)
    plot_w = width - label_w - 40
    y, out, trs = 2, [], []
    for ph, gl in groups.items():
        out.append(f'<circle cx="{label_w-122+5}" cy="{y+11}" r="4.5" fill="{ph_color(ph)}"/>')
        out.append(svg_text(label_w - 110, y + 15, ph if ph in PHYLUM_SLOTS or ph == "(not called)" else ph + " (grey: other phyla)", "t-sec", "start", 'style="font-weight:600"'))
        y += head_h
        for g, n in gl:
            bw = max(plot_w * n / vmax, 8)
            tip = f"{g} ({ph}): {n} inputs"
            out.append(svg_text(label_w - 8, y + 14, short(g, 20), "t-sec", "end"))
            out.append(f'<path class="tip bar" tabindex="0" data-tip="{esc(tip)}" d="M{label_w},{y+3} h{bw-4:.1f} a4,4 0 0 1 4,4 v10 a4,4 0 0 1 -4,4 h-{bw-4:.1f} z" fill="{ph_color(ph)}"/>')
            out.append(f'<rect class="tip hit" tabindex="-1" data-tip="{esc(tip)}" x="{label_w}" y="{y-1}" width="{plot_w+36}" height="{row_h-2}" fill="transparent"/>')
            out.append(svg_text(label_w + bw + 6, y + 14.5, n, "t-pri"))
            trs.append(f"<tr><td>{ph_dot(ph)} {esc(ph)}</td><td>{esc(g)}</td><td>{n}</td></tr>")
            y += row_h
        y += 6
    out.insert(0, f'<line x1="{label_w}" x2="{label_w}" y1="0" y2="{y}" class="axis"/>')
    table = "<table><thead><tr><th>Phylum</th><th>Genus</th><th>Inputs</th></tr></thead><tbody>" + "".join(trs) + "</tbody></table>"
    legend = '<div class="legend">' + "".join(f'<span>{ph_dot(ph)}{esc(ph)}</span>' for ph in groups if ph != "(not called)") + "</div>"
    return f'<svg viewBox="0 0 {width} {y}" class="chart" role="img" aria-label="Genus counts grouped by phylum">{"".join(out)}</svg>', table, legend


def ani_hist():
    vals = [float(r["ani"]) for r in rows if r["ani"] and r["kind"] == "assembly"]
    lo, hi, bw = 0.70, 1.00, 0.01
    nb = int(round((hi - lo) / bw))
    counts = [0] * nb
    for v in vals:
        counts[min(int((v - lo) / bw), nb - 1) if v >= lo else 0] += 1
    W, H, L, R, T, B = 500, 280, 38, 8, 34, 34
    pw, ph = W - L - R, H - T - B
    ticks = nice_ticks(max(counts), 4)
    ymax = ticks[-1]
    out = []
    for t in ticks:
        y = T + ph - ph * t / ymax
        out.append(f'<line x1="{L}" x2="{W-R}" y1="{y:.1f}" y2="{y:.1f}" class="grid"/>')
        out.append(svg_text(L - 8, y + 4, compact(t), "t-mut", "end"))
    slot = pw / nb
    bwid = max(min(slot - 3, 24), 6)
    for i, c in enumerate(counts):
        x = L + i * slot + (slot - bwid) / 2
        a, b = lo + i * bw, lo + (i + 1) * bw
        tip = f"ANI {a:.2f}–{b:.2f}: {c} inputs"
        hgt = ph * c / ymax
        if c:
            out.append(f'<path class="tip bar" tabindex="0" data-tip="{esc(tip)}" d="M{x:.1f},{T+ph} v-{max(hgt-4,0):.1f} a4,4 0 0 1 4,-4 h{bwid-8:.1f} a4,4 0 0 1 4,4 v{max(hgt-4,0):.1f} z" fill="var(--series-1)"/>')
        out.append(f'<rect class="tip hit" data-tip="{esc(tip)}" x="{L+i*slot:.1f}" y="{T}" width="{slot:.1f}" height="{ph}" fill="transparent"/>')
        if c == max(counts):  # direct-label the extreme only
            out.append(svg_text(x + bwid / 2, T + ph - hgt - 6, c, "t-pri", "middle"))
    out.append(f'<line x1="{L}" x2="{W-R}" y1="{T+ph}" y2="{T+ph}" class="axis"/>')
    for v in (0.70, 0.80, 0.90, 1.00):
        x = L + pw * (v - lo) / (hi - lo)
        out.append(svg_text(x, T + ph + 18, f"{v:.2f}", "t-mut", "end" if v >= 1.0 else "middle"))
    # threshold hairlines, directly labelled (solid, not dashed)
    for v, lab in ((0.7765, "genus ≥ 0.777"), (0.9938, "species ≥ 0.994")):
        x = L + pw * (v - lo) / (hi - lo)
        out.append(f'<line x1="{x:.1f}" x2="{x:.1f}" y1="{T-12}" y2="{T+ph}" class="thr"/>')
        out.append(svg_text(x - 4 if v > 0.9 else x + 4, T - 16, lab, "t-sec", "end" if v > 0.9 else "start"))
    out.append(svg_text(L + pw / 2, H - 2, "ANI of the nearest reference (containment, k = 21)", "t-mut", "middle"))
    table = "<table><thead><tr><th>ANI bin</th><th>Inputs</th></tr></thead><tbody>" + "".join(
        f"<tr><td>{lo+i*bw:.2f}–{lo+(i+1)*bw:.2f}</td><td>{c}</td></tr>" for i, c in enumerate(counts) if c) + "</tbody></table>"
    return f'<svg viewBox="0 0 {W} {H}" class="chart" role="img" aria-label="Histogram of nearest-reference ANI">{"".join(out)}</svg>', table


HEAT = ["#cde2fb", "#9ec5f4", "#5598e7", "#2a78d6", "#1c5cab", "#104281"]
HEAT_BINS = [0.90, 0.95, 0.98, 0.99, 0.995]
HEAT_LAB = ["< 0.90", "0.90–0.95", "0.95–0.98", "0.98–0.99", "0.99–0.995", "≥ 0.995"]


def heatmap():
    ids, labs, M = mat["ids"], mat["labels"], mat["ani"]
    n = len(ids)
    cell, L, T = 34, 188, 26
    W, H = L + cell * n + 8, T + cell * n + 62
    out = []
    for i in range(n):
        out.append(svg_text(L - 8, T + i * cell + cell / 2 + 4, f"{i+1}  {short(labs[i][1] or labs[i][0], 24)}", "t-sec", "end"))
        out.append(svg_text(L + i * cell + cell / 2, T - 8, i + 1, "t-mut", "middle"))
        for j in range(n):
            v = M[i][j]
            b = sum(v >= t for t in HEAT_BINS)
            tip = f"{ids[i]} × {ids[j]}: ANI {v:.4f}"
            out.append(f'<rect class="tip" tabindex="0" data-tip="{esc(tip)}" x="{L+j*cell+1}" y="{T+i*cell+1}" width="{cell-2}" height="{cell-2}" rx="3" fill="{HEAT[b]}"/>')
    ly = T + cell * n + 30
    out.append(svg_text(L, ly - 4, "Pairwise ANI", "t-mut"))
    x = L + 82
    for b in range(6):
        out.append(f'<rect x="{x}" y="{ly-14}" width="12" height="12" rx="2" fill="{HEAT[b]}"/>')
        out.append(svg_text(x + 17, ly - 4, HEAT_LAB[b], "t-sec"))
        x += 17 + len(HEAT_LAB[b]) * 6.4 + 16
    head = "".join(f"<th>{i+1}</th>" for i in range(n))
    table = f"<table class='wide'><thead><tr><th></th>{head}</tr></thead><tbody>" + "".join(
        f"<tr><th>{i+1} {esc(ids[i])}</th>" + "".join(f"<td>{M[i][j]:.3f}</td>" for j in range(n)) + "</tr>" for i in range(n)) + "</tbody></table>"
    return f'<svg viewBox="0 0 {W} {H}" class="chart" role="img" aria-label="Pairwise ANI between inputs">{"".join(out)}</svg>', table


def tophit_bars(qid, qhashes):
    hs = hits.get(qid, [])
    if not hs:
        return "", ""
    W, rh, lw = 560, 40, 218
    pw = W - lw - 128
    out, trs = [], []
    for i, h in enumerate(hs):
        c = float(h["shared"]) / max(qhashes, 1)
        ani = c ** (1 / 21)
        y = i * rh + 4
        bw = max(pw * c, 8)
        sp = h["species"] or h["genus"] or "(unnamed)"
        tip = f"{h['hit_id']} · {sp} · {100*c:.1f}% of the query's hashes found (ANI {ani:.4f}) · {h['shared']} shared hashes"
        out.append(svg_text(lw - 10, y + 12, short(sp, 28), "t-pri", "end"))
        out.append(svg_text(lw - 10, y + 27, short(h["hit_id"], 30), "t-mut", "end"))
        out.append(f'<path class="tip bar" tabindex="0" data-tip="{esc(tip)}" d="M{lw},{y+6} h{bw-4:.1f} a4,4 0 0 1 4,4 v10 a4,4 0 0 1 -4,4 h-{bw-4:.1f} z" fill="{"var(--series-1)" if i == 0 else "var(--deemph)"}"/>')
        out.append(f'<rect class="tip hit" data-tip="{esc(tip)}" x="{lw}" y="{y}" width="{pw+110}" height="{rh-4}" fill="transparent"/>')
        out.append(svg_text(lw + bw + 7, y + 21, f"{100*c:.1f}%" + (f"  ·  ANI {ani:.3f}" if i == 0 else ""), "t-pri"))
        trs.append(f"<tr><td>{i+1}</td><td>{esc(h['hit_id'])}</td><td>{esc(sp)}</td><td>{h['shared']}</td><td>{100*c:.1f}%</td><td>{ani:.4f}</td></tr>")
    out.insert(0, f'<line x1="{lw}" x2="{lw}" y1="0" y2="{rh*len(hs)}" class="axis"/>')
    table = "<table><thead><tr><th>#</th><th>Reference</th><th>Species label</th><th>Shared hashes</th><th>Query found</th><th>ANI</th></tr></thead><tbody>" + "".join(trs) + "</tbody></table>"
    return f'<svg viewBox="0 0 {W} {rh*len(hs)+6}" class="chart" role="img" aria-label="Best matching reference genomes">{"".join(out)}</svg>', table


def abund_chart():
    W, H, L, R, T, B = 860, 210, 44, 10, 14, 34
    pw, ph = W - L - R, H - T - B
    keys = [str(i) for i in range(1, 13)]
    vals = [hist["counts"][k] for k in keys] + [hist["over12"]]
    labels = keys + ["13+"]
    ticks = nice_ticks(max(vals), 4)
    ymax = ticks[-1]
    out = []
    for t in ticks:
        y = T + ph - ph * t / ymax
        out.append(f'<line x1="{L}" x2="{W-R}" y1="{y:.1f}" y2="{y:.1f}" class="grid"/>')
        out.append(svg_text(L - 8, y + 4, compact(t), "t-mut", "end"))
    slot = pw / len(vals)
    bwid = min(slot - 6, 24)
    for i, v in enumerate(vals):
        x = L + i * slot + (slot - bwid) / 2
        kept = i + 1 >= hist["min_abund"]
        lab = labels[i]
        tip = f"seen {lab}×: {v} hashes · {'kept' if kept else 'removed (below minimum abundance ' + str(hist['min_abund']) + ')'}"
        hgt = ph * v / ymax
        fill = "var(--series-1)" if kept else "var(--deemph)"
        out.append(f'<path class="tip bar" tabindex="0" data-tip="{esc(tip)}" d="M{x:.1f},{T+ph} v-{max(hgt-4,0):.1f} a4,4 0 0 1 4,-4 h{bwid-8:.1f} a4,4 0 0 1 4,4 v{max(hgt-4,0):.1f} z" fill="{fill}"/>')
        out.append(f'<rect class="tip hit" data-tip="{esc(tip)}" x="{L+i*slot:.1f}" y="{T}" width="{slot:.1f}" height="{ph}" fill="transparent"/>')
        out.append(svg_text(x + bwid / 2, T + ph + 16, lab, "t-mut", "middle"))
    out.append(f'<line x1="{L}" x2="{W-R}" y1="{T+ph}" y2="{T+ph}" class="axis"/>')
    out.append(svg_text(L + pw / 2, H - 2, "times a sketch hash was seen in the reads", "t-mut", "middle"))
    out.append(svg_text(L + slot * 1.15, T + 14, f"← {compact(vals[0])} hashes seen once: mostly sequencing errors, removed", "t-sec", "start"))
    table = "<table><thead><tr><th>Seen</th><th>Hashes</th><th>Used</th></tr></thead><tbody>" + "".join(
        f"<tr><td>{labels[i]}×</td><td>{vals[i]}</td><td>{'yes' if i+1 >= hist['min_abund'] else 'no'}</td></tr>" for i in range(len(vals))) + "</tbody></table>"
    return f'<svg viewBox="0 0 {W} {H}" class="chart" role="img" aria-label="Abundance histogram of read hashes">{"".join(out)}</svg>', table


# ------------------------------------------------------------------ components
ICON = {  # icon + label always travel together; colour is never the only channel
    "ok": '<svg class="ic" viewBox="0 0 16 16" aria-hidden="true"><circle cx="8" cy="8" r="7" fill="var(--st-good)"/><path d="M4.6 8.3l2.2 2.2 4.5-4.8" fill="none" stroke="#fff" stroke-width="1.8" stroke-linecap="round" stroke-linejoin="round"/></svg>',
    "low_confidence": '<svg class="ic" viewBox="0 0 16 16" aria-hidden="true"><path d="M8 1.5l7 12.5H1z" fill="var(--st-warn)"/><path d="M8 6v3.6M8 11.6v.4" stroke="#0b0b0b" stroke-width="1.7" stroke-linecap="round"/></svg>',
    "no_match": '<svg class="ic" viewBox="0 0 16 16" aria-hidden="true"><circle cx="8" cy="8" r="7" fill="var(--st-serious)"/><path d="M5 11L11 5" stroke="#0b0b0b" stroke-width="1.8" stroke-linecap="round"/></svg>',
    "error": '<svg class="ic" viewBox="0 0 16 16" aria-hidden="true"><circle cx="8" cy="8" r="7" fill="var(--st-crit)"/><path d="M5.3 5.3l5.4 5.4M10.7 5.3l-5.4 5.4" stroke="#fff" stroke-width="1.8" stroke-linecap="round"/></svg>',
}
STATUS_LABEL = {"ok": "classified", "low_confidence": "low confidence", "no_match": "no match", "error": "error"}


def chip(st):
    return f'<span class="chip">{ICON[st]}{STATUS_LABEL[st]}</span>'


def reasons(r):
    """Why was each rank withheld? Mirrors the call logic with the 'known' thresholds."""
    out = []
    best = RANKS.index(r["best_rank"]) if r["best_rank"] else -1
    ani = float(r["ani"]) if r["ani"] else 0.0
    for i, k in enumerate(RANKS):
        if i <= best:
            continue
        cons = float(r[f"conf_{k}"] or 0)
        if not r[f"call_{k}"] and (cons == 0 or not r["top_hit"]):
            out.append((k, f"The nearest reference genomes have no {k} recorded in the reference metadata, so there is nothing to call. This reflects the reference annotation, not weak evidence"))
        elif cons < TH["min_cons"][i]:
            out.append((k, f"Vote consistency {cons:.3f} is below the {TH['min_cons'][i]:.3f} required: the near-best references disagree"))
        elif ani < TH["min_ani"][i]:
            out.append((k, f"ANI {ani:.4f} is below the {TH['min_ani'][i]:.4f} required for a {k} call"))
        else:
            out.append((k, "a higher rank was withheld, so deeper ranks are not called"))
        break  # deeper ranks follow from the first withheld one
    return out


def ladder(r):
    best = RANKS.index(r["best_rank"]) if r["best_rank"] else -1
    cells = []
    for i, k in enumerate(RANKS):
        called = i <= best and r[f"call_{k}"]
        cons = float(r[f"conf_{k}"] or 0)
        gap = (not called) and i < best  # a deeper rank was called, this one has no value in the reference metadata
        name = r[f"call_{k}"] if called else ("not in reference" if gap else "withheld")
        tip = (f"{k}: {r[f'call_{k}']} · vote consistency {cons:.2f}" if called else (f"{k}: the reference metadata has no {k} for these genomes" if gap else f"{k}: not called"))
        cells.append(
            f'<div class="rk {"on" if called else ("na" if gap else "off")} tip" tabindex="0" data-tip="{esc(tip)}"{(" style=\"border-left-color:" + ph_color(r["call_phylum"]) + "\"") if (k == "phylum" and called) else ""}><div class="rk-n">{k}</div>'
            f'<div class="rk-v">{esc(name)}</div><div class="rk-m"><span style="width:{100*cons if called else 0:.0f}%"></span></div></div>')
    return '<div class="ladder">' + "".join(cells) + "</div>"


def stat(label, value, st=None):
    ic = ICON[st] if st else ""
    return f'<div class="tile"><div class="tile-l">{ic}{esc(label)}</div><div class="tile-v">{value}</div></div>'


def qc_panel(r):
    tb = float(r["total_bases"] or 0)
    items = [("Bases", f"{tb/1e6:.1f} Mb" if tb >= 1e6 else compact(tb)), ("N50", compact(float(r["n50"])) if r["n50"] not in ("", "0") else "n/a (reads)"),
             ("Sketch hashes", r["n_hashes"])]
    flags = [x for x in r["qc_flags"].split(",") if x]
    fl = "".join(f'<span class="flag">{esc(x)}</span>' for x in flags) or '<span class="muted">none</span>'
    return '<div class="qc">' + "".join(f'<div><span class="muted">{a}</span><b>{esc(str(b))}</b></div>' for a, b in items) + f'<div><span class="muted">QC flags</span>{fl}</div></div>'


def card(r, why):
    qid = r["query_id"]
    st = r["status"]
    head = f'<summary><span class="card-id">{esc(qid)}</span>{chip(st)}{ph_dot(r["call_phylum"]) if r["call_phylum"] else ""}<span class="card-call">{esc(r["call_species"] or r["call_genus"] or "")}</span></summary>'
    if st == "error":
        return (f'<details class="card" id="q-{esc(qid)}" open>{head}<div class="card-b"><p class="why">This input could not be processed: '
                f'<b>{esc(r["message"])}</b>. The run continued; this row stays in results.tsv with status <code>error</code>.</p></div></details>')
    qh = int(r["n_hashes"] or 1)
    bars, bars_t = tophit_bars(qid, qh)
    ev = "".join(f"<li><b>{k.capitalize()} not called.</b> {esc(t)}.</li>" for k, t in reasons(r))
    if st == "no_match":
        ev = "<li><b>No rank called.</b> The nearest reference shares too little with this input to support even a phylum call. Not a fungus, or a lineage the reference does not cover.</li>"
    ev_html = f'<ul class="why">{ev}</ul>' if ev else '<p class="why ok">Every rank down to species passed its threshold.</p>'
    clus = f'<div><span class="muted">ANI cluster of the top hit</span><b class="mono">{esc(short(r["cluster"], 40)) or "—"}</b></div>' if r["cluster"] else ""
    reads = ""
    if r["kind"] == "reads" and qid.startswith("A1163_5x"):
        ch, tb = abund_chart()
        reads = (f'<h4>Reads: what was kept</h4><figure class="fig">{figtools(ch, tb)}<figcaption>Hashes seen once are mostly sequencing errors and are dropped '
                 f'(minimum abundance {hist["min_abund"]}). The peak of the kept hashes (median {hist["median_kept_abundance"]}×) estimates k-mer coverage; '
                 f'{100*float(r["containment_query"]):.1f}% of the kept hashes are found in the best reference.</figcaption></figure>')
    return (f'<details class="card" id="q-{esc(qid)}" {"open" if qid in OPEN else ""}>{head}<div class="card-b">'
            f'<h4>Call by rank</h4>{ladder(r)}{ev_html}'
            f'<div class="two"><div><h4>Best matching references</h4><figure class="fig">{figtools(bars, bars_t)}<figcaption>Bar = share of this input\'s sketch hashes found in the reference; '
            f'the first bar is the top hit.</figcaption></figure></div><div><h4>Input</h4>{qc_panel(r)}<div class="qc">'
            f'<div><span class="muted">ANI to top hit</span><b>{f(r["ani"])}</b></div><div><span class="muted">Shared hashes</span><b>{r["shared"]}</b></div>{clus}</div></div></div>{reads}</div></details>')


def figtools(svg, table):
    return f'<div class="figbtn"><button class="tbtn" data-table-toggle aria-pressed="false">Table view</button></div><div class="chartbox">{svg}</div><div class="tableview" hidden>{table}</div>'


# ------------------------------------------------------------------ choose example cards (all real inputs)
def pick(pred, key=lambda r: 0, rev=False):
    c = [r for r in rows if pred(r)]
    return sorted(c, key=key, reverse=rev)[0] if c else None


ex_species = pick(lambda r: r["status"] == "ok" and r["best_rank"] == "species" and r["kind"] == "assembly" and "Aspergillus" in r["call_species"], lambda r: -float(r["ani"] or 0))
ex_genus = pick(lambda r: r["status"] == "ok" and r["best_rank"] == "genus" and r["kind"] == "assembly" and float(r["ani"] or 0) > 0.95, lambda r: -float(r["ani"]))
ex_low = pick(lambda r: r["status"] == "low_confidence")
ex_none = pick(lambda r: r["status"] == "no_match")
ex_reads = byid["A1163_5x"]
ex_err = pick(lambda r: r["status"] == "error")
EX = [x for x in (ex_species, ex_genus, ex_low, ex_none, ex_reads, ex_err) if x]
OPEN = {ex_species["query_id"], ex_reads["query_id"]}

# ------------------------------------------------------------------ assemble the page
depth_svg, depth_tbl = depth_chart()
ph_ct = collections.Counter(r["call_phylum"] or "(not called)" for r in rows if r["status"] != "error")
ge_ct = collections.Counter(r["call_genus"] for r in rows if r["call_genus"])
ph_items = sorted(ph_ct.items(), key=lambda kv: (-kv[1], kv[0]))
ph_svg, ph_tbl = hbars(ph_items, "Phylum", unit=" inputs", colors=[ph_color(n) for n, _ in ph_items])
genus_ph = {}
for r in rows:
    if r["call_genus"] and r["status"] != "error":
        genus_ph.setdefault(r["call_genus"], collections.Counter())[r["call_phylum"] or "(not called)"] += 1
ge_items = [(g, genus_ph[g].most_common(1)[0][0], n) for g, n in ge_ct.most_common(12)]
ge_svg, ge_tbl, ge_legend = genus_by_phylum(ge_items)
ani_svg, ani_tbl = ani_hist()
heat_svg, heat_tbl = heatmap()
pct = round(100 * classified / ok_total)
_low = sorted([r for r in rows if r["ani"] and r["status"] != "error"], key=lambda r: float(r["ani"]))[:5]
low_tbl = '<table class="small full"><thead><tr><th>Input</th><th>Call</th><th class="num">ANI</th></tr></thead><tbody>' + "".join(
    f'<tr><td class="mono">{esc(short(r["query_id"], 22))}</td><td>{esc(r["call_species"] or r["call_genus"] or r["call_family"] or r["call_phylum"] or "no call")}</td><td class="num">{float(r["ani"]):.3f}</td></tr>' for r in _low) + "</tbody></table>"

idx = run["index"]
thr_rows = "".join(
    f"<tr><td>{k}</td><td>{TH['min_ani'][i]:.4f}</td><td>{TH['min_cons'][i]:.2f}</td></tr>" for i, k in enumerate(RANKS))


def trow(r):
    st = r["status"]
    call = r["call_species"] or r["call_genus"] or r["call_family"] or r["call_order"] or r["call_class"] or r["call_phylum"] or ""
    has_card = r["query_id"] in {x["query_id"] for x in EX}
    link = f'<a href="#q-{esc(r["query_id"])}">{esc(r["query_id"])}</a>' if has_card else esc(r["query_id"])
    best = r["best_rank"] or "—"
    flags = "".join(f'<span class="flag">{esc(x)}</span>' for x in r["qc_flags"].split(",") if x)
    return (f'<tr data-status="{st}" data-q="{esc((r["query_id"]+" "+call+" "+r["call_phylum"]).lower())}"><td data-v="{esc(r["query_id"])}">{link}</td><td data-v="{st}">{chip(st)}</td>'
            f'<td data-v="{esc(r["call_phylum"])}">{ph_dot(r["call_phylum"]) if r["call_phylum"] else ""}{esc(r["call_phylum"]) or "—"}</td>'
            f'<td data-v="{RANKS.index(r["best_rank"]) if r["best_rank"] else -1}">{best}</td><td data-v="{esc(call)}">{esc(call) or "—"}</td>'
            f'<td data-v="{r["ani"] or 0}" class="num">{f(r["ani"],3) or "—"}</td><td>{flags}</td></tr>')


table_html = "".join(trow(r) for r in rows)
cards_html = "".join(card(r, None) for r in EX)

CSS = open(os.path.join(HERE, "mockup.css")).read()
JS = open(os.path.join(HERE, "mockup.js")).read()

page = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>GmG report (mockup)</title><style>{CSS}</style></head>
<body>
<div class="mockbar"><b>DESIGN MOCKUP (issue #22).</b> Real output of <code>gmg run</code> on 148 inputs; layout and wording are what is being reviewed. Charts: hover for values, use “Table view” for the numbers.</div>
<main class="page viz-root">
<header class="top">
  <div><div class="eyebrow">GmG · Guess My Genome</div><h1>Classification report</h1>
  <p class="sub">Run {esc(run['started_utc'][:10])} · {n_in} inputs · gmg {esc(run['version'])} ({esc(run['commit'][:9])}) · index of {idx['n_genomes']:,} genomes (k = {idx['k']}, scaled = {idx['scaled']}) · thresholds: {esc(run['thresholds']['name'])}</p></div>
  <button id="theme" class="tbtn" aria-label="Toggle dark mode">Dark / light</button>
</header>

<section class="summary" aria-label="Summary">
  <div class="hero"><div class="hero-v">{pct}%</div><div class="hero-l">of the {ok_total} processed inputs were classified to genus or species</div></div>
  <div class="tiles">
    {stat('species', depth_ct['species'])}{stat('genus', depth_ct['genus'])}
    {stat('low confidence', depth_ct['family'] + depth_ct['higher'], 'low_confidence')}
    {stat('no match', depth_ct['none'], 'no_match')}{stat('errors', status_ct['error'], 'error')}
  </div>
</section>

<section class="card2"><h2>How deep each call went</h2>
<p class="lede">A rank is called only when its evidence passes the calibrated threshold; otherwise the report stops higher up and says why.</p>
<figure class="fig">{figtools(depth_svg, depth_tbl)}</figure></section>

<section class="grid2">
  <div class="card2"><h2>What the run contains</h2><h3>Phylum</h3><figure class="fig">{figtools(ph_svg, ph_tbl)}<figcaption>Colours are fixed per phylum (the eight largest phyla in the reference). Other phyla share grey; the table view names every phylum.</figcaption></figure>
  <h3>Genus (top {len(ge_items)}), coloured by phylum</h3><figure class="fig">{figtools(ge_svg, ge_tbl)}{ge_legend}</figure></div>
  <div class="card2"><h2>How close the nearest reference is</h2>
  <p class="lede">Almost every input has a near-identical reference. Calls to species need ANI ≥ {TH['min_ani'][5]:.3f}; the few lower values are less well covered lineages.</p>
  <figure class="fig">{figtools(ani_svg, ani_tbl)}</figure>
  <h3>Lowest ANI to a reference: look at these first</h3>{low_tbl}</div>
</section>

<section class="card2"><h2>Relationships among inputs <span class="opt">optional · up to 60 inputs</span></h2>
<p class="lede">Pairwise ANI between 18 of the inputs, computed from the same sketches. Blocks along the diagonal are species and species complexes.</p>
<figure class="fig">{figtools(heat_svg, heat_tbl)}</figure></section>

<section class="card2" id="results"><h2>All inputs</h2>
<div class="filters"><label>Status <select id="f-status"><option value="">all</option><option value="ok">classified</option><option value="low_confidence">low confidence</option><option value="no_match">no match</option><option value="error">error</option></select></label>
<label>Search <input id="f-q" type="search" placeholder="input or taxon"></label><span id="f-n" class="muted"></span>
<a class="dl" href="sample/results.tsv" download>Download results.tsv</a></div>
<div class="tscroll"><table class="results"><thead><tr><th data-sort>Input</th><th data-sort>Status</th><th data-sort>Phylum</th><th data-sort>Deepest rank</th><th data-sort>Call</th><th data-sort class="num">ANI</th><th>QC flags</th></tr></thead>
<tbody>{table_html}</tbody></table></div>
<p class="muted small">Only inputs with a detail card below are linked in this mockup; the real report has a card for every input, collapsed.</p></section>

<section><h2 class="cards-h">Per-input detail</h2>{cards_html}</section>

<section class="card2 meth"><h2>How to read this report</h2>
<ul>
<li><b>Call by rank.</b> Phylum to species. A rank is called when the nearest reference genomes agree on it (vote consistency) and the input is similar enough to them (ANI). Higher ranks are easier than lower ones, so the ladder can stop part-way.</li>
<li><b>Vote consistency</b> is the share of near-best references agreeing on the taxon, an evidence value, not a probability.</li>
<li><b>ANI</b> here is containment ANI: the k = 21 sketch fraction of the input found in a reference, raised to 1/21.</li>
<li><b>Species labels follow the reference metadata.</b> Where references disagree about a name (hybrids, species complexes) the genus is usually right and the species is withheld; the ANI cluster is the safer unit.</li>
<li>Thresholds are calibrated to about 95% precision on held-out genomes of the reference (<code>known</code> set). Genomes from a genus the reference lacks need the stricter <code>novel</code> set.</li>
</ul>
<h3>Thresholds used ({esc(run['thresholds']['name'])})</h3>
<table class="small"><thead><tr><th>Rank</th><th>Minimum ANI</th><th>Minimum vote consistency</th></tr></thead><tbody>{thr_rows}</tbody></table>
<h3>Reproducibility</h3>
<p class="small mono">{esc(' '.join(run['command'][:1]))} run … · gmg {esc(run['version'])} ({esc(run['commit'])}) · index {esc(os.path.basename(idx['path']))} ({idx['bytes']:,} bytes) · {run['elapsed_seconds']} s · results.tsv · hits.tsv · run.json</p>
</section>
</main>
<div id="tooltip" role="tooltip" hidden></div>
<script>{JS}</script>
</body></html>"""
open(os.path.join(HERE, "mockup.html"), "w").write(page)
print("wrote mockup.html", len(page) // 1024, "KB;", "cards:", [x["query_id"] for x in EX])
