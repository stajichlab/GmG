#!/usr/bin/env python3
"""Consolidate every genome / read set whose label looks wrong into one hand-off table (#14).

usage: make_upstream_report.py [OUT_DIR=reports/label_audit]
Writes OUT_DIR/upstream_report.tsv (machine readable) and the tables of OUT_DIR/upstream_report.md (between the markers
<!-- TABLES --> ... <!-- /TABLES -->; the prose around them is hand-written and kept).

Inputs (all in the repository): reports/label_audit/{reference_exclusions.txt,exclude_ids.txt,loo_suspects.tsv,hub_references.tsv,
size_outliers.tsv,log_conflicts.tsv}, data/external/bfd_meta.tsv (BFD label of every reference genome), data/validation/
{truth_set.list.tsv,disagreements.tsv,sra_runs.tsv}. The truth-set disagreements of the FIRST run (index 2026.10.2, commit 374dfec) are
read from git, because the genomes of the 14 later exclusions are no longer neighbours in the re-run. The rows about lab/provider
labels (ZyGoLife, IndigoAg, Kasson reads and assemblies) come from per-genome results that stay on the HPCC (GmG_runs/external*,
GmG_runs/reads/titration, GmG_runs/reads_sra) and are typed in EXTERNAL below with the file they were read from.

Tier rules (also in the report): A = at least two independent lines of evidence and ANI >= 0.98 (assemblies or reads of the same
strain agree, or hub + leave-one-out, or hub/leave-one-out + a RefSeq genome); B = one line with ANI >= 0.98, or two lines that do not
reach A; C = weak or ambiguous (uninformative 'sp.' label, order-level only, ANI < 0.98, size outlier only, metagenome input, taxonomy
update). Suspicions come from sequence similarity only.
"""
import csv, collections, subprocess, sys, os, re

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, 'reports/label_audit')


def rd(p):
    return list(csv.DictReader(open(os.path.join(ROOT, p) if not p.startswith('/') else p), delimiter='\t'))


meta = {r['id']: r for r in rd('data/external/bfd_meta.tsv')}
ex57 = [l.strip() for l in open(os.path.join(ROOT, 'reports/label_audit/reference_exclusions.txt')) if l.strip()]
ex43 = set(open(os.path.join(ROOT, 'reports/label_audit/exclude_ids.txt')).read().split())
EX = set(ex57)
loo = rd('reports/label_audit/loo_suspects.tsv')
hubs = {r['top_id']: r for r in rd('reports/label_audit/hub_references.tsv')}
sizeo = rd('reports/label_audit/size_outliers.tsv')
logc = rd('reports/label_audit/log_conflicts.tsv')
truth = {r['id']: r for r in rd('data/validation/truth_set.list.tsv')}
first = subprocess.run(['git', '-C', ROOT, 'show', '374dfec:data/validation/disagreements.tsv'], capture_output=True, text=True, check=True).stdout
dis = list(csv.DictReader(first.splitlines(), delimiter='\t'))
sra = {r['run']: r for r in rd('data/validation/sra_runs.tsv')}

RANKS = ['phylum', 'class', 'order', 'family', 'genus']
acc = lambda i: '_'.join(i.split('_')[:2])
lin = lambda m: ';'.join(m[k] for k in RANKS) if m else ''
label = lambda i: meta[i]['species'] if i in meta else ''

COLS = ['tier', 'row_id', 'kind', 'accession', 'bfd_id', 'current_organism', 'current_lineage', 'best_hit_accession', 'best_hit_taxon', 'ani',
        'size_ratio', 'flags', 'evidence_type', 'evidence', 'suggested_taxon', 'suggested_action', 'in_reference_2026.10.3', 'sources']
rows = []


def add(**k):
    r = {c: '' for c in COLS}
    r.update(k)
    rows.append(r)


# ---------------------------------------------------------------- 57 genomes left out of the reference
import gzip
clus = {r['id']: r for r in csv.DictReader(gzip.open(os.path.join(ROOT, 'data/clusters/bfd_ani0.95.tsv.gz'), 'rt'), delimiter='\t')}
cmembers = collections.defaultdict(list)
for _i, _r in clus.items():
    cmembers[_r['cluster']].append(_i)


def truth_lines(g):
    """RefSeq genomes (first run) whose nearest hit was g at ANI >= 0.95 with a different lineage."""
    d = collections.OrderedDict()
    for x in dis:
        if x['top_hit'] == g and x['category'] == 'reference_label_error' and float(x['ani'] or 0) >= 0.95:
            d.setdefault(x['id'], []).append(x)
    return d


meta_class = {}
for _m in meta.values():
    meta_class.setdefault(_m['species'], _m['class'])


def plural(counter, n):
    """label of the plurality species (or its genus) among n others"""
    if not counter:
        return ''
    sp, c = counter.most_common(1)[0]
    if c * 2 >= n and ' sp.' not in sp:
        return f"{sp} ({c} of {n})"
    gc = collections.Counter(s.split()[0] for s in counter.elements())
    gname, gcount = gc.most_common(1)[0]
    if gcount * 2 >= n and gname[0].isupper() and gname.lower() not in ('uncultured',):
        return f"{gname} ({gcount} of {n}; e.g. {sp})"
    cc = collections.Counter(meta_class.get(s, '') for s in counter.elements())
    cname, ccount = cc.most_common(1)[0]
    return f"{cname} ({ccount} of {n}; e.g. {sp})"


def assess(g, kept=False):
    m = meta[g]
    own = [x for x in (loo + logc if kept else loo) if x['id'] == g and (not kept or x['top_id'] not in EX)]
    pointing = [x for x in loo + logc if x['top_id'] == g and (not kept or x['id'] not in EX)]
    ph = [x for x in own if x['rank'] == 'phylum']
    cl_ = [x for x in own if x['rank'] == 'class']
    best = (ph or cl_ or own or [None])[0]
    h = hubs.get(g)
    if kept and not h:
        pid = set(x['id'] for x in pointing)
        if len(pid) >= 3:
            h = {'n': str(len(pid)), 'their_label': collections.Counter(x['own'] for x in pointing).most_common(1)[0][0], 'neighbours_say': pointing[0]['neighbours']}
    tl = truth_lines(g)
    ani_own = max([float(x['ani']) for x in own] or [0])
    ani_pt = max([float(x['ani']) for x in pointing] or [0])
    tani = max([float(x['ani']) for v in tl.values() for x in v] or [0])
    ani_all = max(ani_own, ani_pt, tani)
    # cluster at 95% ANI: do the other members share this genome's label?
    if kept:
        rank = 'class' if (any(x['rank'] in ('phylum', 'class') for x in own) or (h and any(x['rank'] in ('phylum', 'class') for x in pointing))) else 'order'
    else:
        rank = 'class' if (own or h) else 'genus'
    others = [o for o in cmembers[clus[g]['cluster']] if o != g]
    same = sum(1 for o in others if meta[o][rank] == m[rank] and m[rank])
    ncl = len(others)
    _oc = collections.Counter(meta[o][rank] for o in others if meta[o][rank] != m[rank])
    coherent = bool(_oc) and _oc.most_common(1)[0][1] >= 2 and _oc.most_common(1)[0][1] * 2 >= ncl
    odd = ncl >= 2 and same <= 0.05 * ncl and bool(m[rank]) and coherent
    supported = same >= 2 and same * 2 >= ncl
    # the main neighbour's own cluster: if >= 2 of its cluster mates share ITS label, the neighbour is not the odd one
    nb0 = own[0]['top_id'] if own else None
    odd_pt = False
    if nb0 and nb0 in clus and not odd and not supported:
        pcl = [o for o in cmembers[clus[nb0]['cluster']] if o != nb0 and o != g]
        psame = sum(1 for o in pcl if meta[o][rank] == meta[nb0][rank] and meta[nb0][rank])
        odd_pt = len(pcl) >= 2 and psame >= 2 and psame * 2 >= len(pcl) and bool(meta[nb0]['genus']) and meta[nb0][rank] != m[rank]
    ocount = collections.Counter(meta[o]['species'] for o in others if meta[o][rank] != m[rank])
    ev, et = [], []
    if own:
        et.append('LOO')
        nb = own[0]['top_id']
        ev.append(f"leave-one-out: nearest genome {acc(nb)} ({label(nb)}) at ANI {ani_own:.4f} has a different {'/'.join(dict.fromkeys(x['rank'] for x in own))} ({best['own']} here, {best['neighbours']} there)"
                  + ('; that genome is also left out' if nb in EX else ''))
    if h:
        et.append('HUB')
        ev.append(f"hub: nearest genome of {h['n']} other BFD genomes labelled {h['their_label']} (this one: {h['neighbours_say']})")
    if tl:
        et.append('REFSEQ')
        names = [f"{k} {truth[k]['organism'] if k in truth else ''} (ANI {max(float(x['ani']) for x in v):.4f})" for k, v in list(tl.items())[:3]]
        ev.append(f"{len(tl)} RefSeq genome(s) of another {'genus' if rank == 'genus' else 'lineage'} have it as nearest neighbour at ANI >= 0.95: " + '; '.join(names) + (' ...' if len(tl) > 3 else ''))
    if ncl:
        et.append('CLUSTER')
        if odd:
            ev.append(f"odd one out of its 95%-ANI cluster: none of the {ncl} other genomes shares its {rank} (they are {plural(ocount, ncl)})")
        elif supported:
            ev.append(f"its {rank} is shared by {same} of the {ncl} other genomes of its 95%-ANI cluster: the label is supported, the conflict is on the neighbour's side")
        else:
            ev.append(f"95%-ANI cluster of {ncl + 1}: {same} of {ncl} others share its {rank}; side undetermined")
    else:
        ev.append('no other genome in its 95%-ANI cluster, so which side of the pair is wrong is not determined')
    # suggested taxon
    if odd and ocount:
        sug = plural(ocount, ncl)
    elif tl:
        gen = collections.Counter(truth[k]['genus'] for k in tl if k in truth)
        sug = gen.most_common(1)[0][0] if gen else ''
    elif pointing:
        sp = collections.Counter(label(x['id']) for x in pointing)
        sug = f"{h['their_label'] if h else pointing[0]['own']} (e.g. {sp.most_common(1)[0][0]})"
    elif best:
        sug = label(best['top_id'])
    else:
        sug = ''
    uninformative = (m['genus'] == '')
    partner_hub = bool(own) and all(x['top_id'] in hubs for x in own) and not h
    lines = (1 if odd else 0) + (1 if h and int(h['n']) >= 3 else 0) + (1 if tl else 0) + (1 if own else 0)
    if uninformative:
        tier = 'C'
    elif supported:
        tier = 'C'
    elif (odd or odd_pt) and (h or tl or ani_all >= 0.98):
        tier = 'A'
    elif ani_all >= 0.98 or lines >= 2:
        tier = 'B'
    else:
        tier = 'C'
    if ph and tier == 'C' and not uninformative and not supported:
        tier = 'B'  # a phylum-level conflict at ANI >= 0.95 is never only a taxonomy difference
    if partner_hub and not odd and tier == 'B':
        tier = 'C'
    # a specific name whose only conflicting neighbours carry no genus is probably the better-labelled side of the pair
    if own and all(not meta.get(x['top_id'], {}).get('genus', 'x') for x in own) and not odd and not h:
        tier = 'C'
        supported = True
    if odd_pt and tier == 'A':
        ev[:] = [e for e in ev if 'side undetermined' not in e and 'which side' not in e]
        ev.append(f"its nearest genome {acc(nb0)} ({label(nb0)}) is itself in a 95%-ANI cluster whose other members share its {rank}, so the neighbour is not the odd one")
        if not odd:
            sug = label(nb0)
    nbk = [x['top_id'] for x in own if x['top_id'] not in EX]
    if tier == 'A':
        act = f"Relabel (or flag as misidentified/contaminated): sequence is that of {sug}; check source material and strain; if unresolved, remove."
    elif uninformative:
        act = f"Keep, optionally refine: the label is only '{m['species']}', the nearest genomes name {sug}. Not a wrong name."
    elif supported:
        act = "Keep: probably the better-labelled side (its label is shared by the genomes of its cluster, or the conflicting neighbour has no genus in its label)."
    elif tier == 'B':
        act = f"Check: this genome or its neighbour carries a wrong name or is contaminated; similarity cannot say which. If this one: relabel to {sug}."
    else:
        act = f"Keep unless checked: ANI below 0.98 or a single weak line; possible taxonomy difference. Candidate: {sug}."
    flags = ['kept_in_reference'] if kept else ['audit43' if g in ex43 else 'truthset14']
    if h:
        flags.append(f"hub_n={h['n']}")
    flags.append(f"cluster95_others={ncl}")
    if nbk:
        flags.append('neighbour_in_reference')
    bacc = acc(best['top_id']) if best else (acc(next(iter(tl))) if tl else '')
    btax = label(best['top_id']) if best else (truth[next(iter(tl))]['organism'] if tl and next(iter(tl)) in truth else '')
    add(tier=tier, kind='bfd_genome', accession=acc(g), bfd_id=g, current_organism=m['species'], current_lineage=lin(m), best_hit_accession=bacc, best_hit_taxon=btax,
        ani=f"{ani_all:.4f}" if ani_all else '', flags=';'.join(flags), evidence_type='+'.join(et), evidence=' | '.join(ev),
        suggested_taxon=sug, suggested_action=act, **{'in_reference_2026.10.3': 'yes' if kept else 'no (left out)'},
        sources=('reports/label_audit/{loo_suspects,log_conflicts}.tsv; data/clusters/bfd_ani0.95.tsv.gz' if kept else 'reports/label_audit/{exclude_ids,loo_suspects,hub_references}; reference_exclusions.txt; data/clusters/bfd_ani0.95.tsv.gz' + ('' if g in ex43 else '; first truth-set run (disagreements.tsv @374dfec); truth_set_reviews.tsv')))


for g in ex57:
    assess(g)

# ---------------------------------------------------------------- genomes that stay in the reference but conflict with a genome that stays too
# (a genome whose conflicting neighbour is one of the 57 is a victim of that neighbour and not listed)
kept_ids = collections.OrderedDict()
for x in loo + logc:
    if x['id'] not in EX and x['top_id'] not in EX:
        kept_ids[x['id']] = 1
_pt = collections.defaultdict(set)
for x in loo + logc:
    if x['id'] not in EX and x['top_id'] not in EX:
        _pt[x['top_id']].add(x['id'])
for t, v in _pt.items():
    if len(v) >= 3:
        kept_ids[t] = 1
for g in kept_ids:
    if g in meta:
        assess(g, kept=True)

# ---------------------------------------------------------------- size outliers (none is left out)
for r in sizeo:
    m = meta.get(r['id'])
    add(tier='C', kind='bfd_genome', accession=acc(r['id']), bfd_id=r['id'], current_organism=r['sp'], current_lineage=lin(m), size_ratio=r['x_median'],
        flags='size_outlier', evidence_type='SIZE', evidence=f"sketch size {r['x_median']}x the median of its species (species with >= 5 genomes; {r['n_hashes']} hashes): possible mixed or contaminated assembly, or a true polyploid or hybrid",
        suggested_taxon='', suggested_action='Keep; submitter could check the assembly for contamination or duplicated haplotypes.', **{'in_reference_2026.10.3': 'yes'},
        sources='reports/label_audit/size_outliers.tsv')

# ---------------------------------------------------------------- the SRA run
add(tier='B', kind='sra_run+bfd_genome', accession='ERR14883901 / SAMEA118140668 / GCA_982397225.1', bfd_id='GCA_982397225.1_27-FMI', current_organism='Coccidioides immitis',
    current_lineage=lin(meta['GCA_982397225.1_27-FMI']), best_hit_accession='GCA_982397405.1', best_hit_taxon='Ophidiomyces ophidiicola (T55-F)', ani='0.9473 (reads, own genome hidden); 0.9503 (assembly, leave-one-out)',
    flags='5 other C. immitis genomes in the reference; next hit 19 shared hashes vs 772', evidence_type='SRA+LOO',
    evidence='reads of BioSample SAMEA118140668 (own genome hidden) and the assembly GCA_982397225.1 both have Ophidiomyces ophidiicola GCA_982397405.1 as best hit (772 shared hashes, next hit 19), not any of the 5 other C. immitis; the pair is reciprocal (T55-F points back at 27-FMI at 0.9503). Same BioSample, so the two lines are not fully independent. Neither genome is among the 57 left out, both are in the reference.',
    suggested_taxon='not C. immitis; Onygenaceae, nearest Ophidiomyces ophidiicola (ANI 0.95 is below species level, so the genus is not certain)',
    suggested_action='Check: this genome or GCA_982397405.1 (T55-F, same project, labelled Ophidiomyces ophidiicola) carries a wrong name; side undetermined.',
    **{'in_reference_2026.10.3': 'yes'}, sources='reports/validation/reads_sra.md; data/validation/sra_runs.tsv; GmG_runs/reads_sra/*_n500000/results.tsv; GmG_runs/bench/bench_loo_s10000.tsv')

# ---------------------------------------------------------------- decoy
add(tier='C', kind='decoy_assembly', accession='GCF_965202585.1', bfd_id='(decoy, not in BFD)', current_organism='Thalia democratica (salp)', current_lineage='Chordata',
    best_hit_accession='9 ZyGoLife assemblies + 1 read set', best_hit_taxon='Conidiobolus/Neoconidiobolus (7), Mucor hiemalis (2), Conidiobolus antarcticus ARSEF 691 reads',
    ani='0.83-0.913 (best hit of 9 assemblies, 1 read set)', flags='false not_fungal on Entomophthoromycota (2 assemblies at ANI 0.9001/0.9055, 1 read set at 0.913)', evidence_type='ASM+READS',
    evidence='this animal genome shares 825 and 982 hashes (11-12% of the assembly) with two Conidiobolus assemblies and 182 (15%) with one read set, closer (ANI 0.90-0.91) than any fungal genome in the reference (0.88-0.90); chance cannot give a salp 10% of the k-mers of a fungus. Not examined: whether the assembly holds fungal contigs. Hiding it removed all three false calls (scripts/hpcc/28_thalia_check.sbatch).',
    suggested_taxon='unknown; possible fungal contaminant contigs in the assembly, or a conserved repeat', suggested_action='Check (submitter or NCBI contamination screen, e.g. FCS-GX): does GCF_965202585.1 contain Entomophthoromycota contigs? Removed from the GmG decoy panel in 2026.10.4.',
    **{'in_reference_2026.10.3': 'n/a (decoy; removed in 2026.10.4)'}, sources='reports/validation/decoys.md (addendum); reports/validation/reads_titration.md; reports/validation/external_zygolife.md')

# ---------------------------------------------------------------- lab / provider labels (not NCBI accessions)
LAB = [
    # tier, kind, sample, label, best_hit, taxon, ani, flags, evidence_type, evidence, suggested, action, source
    ('A', 'lab_read_set+lab_assembly', 'NRRL A-19777', 'Absidia sp. NRRL A-19777', 'GCA_000945115.1', 'Lichtheimia ramosa', '0.9991 reads; 1.0000 assembly', 'ok (both)', 'READS+ASM',
     'reads (1.5M pairs) and assembly of the same strain both have Lichtheimia ramosa as best hit (2700 shared hashes in the reads)', 'Lichtheimia ramosa', 'Relabel in the lab sheet (Absidia sp. -> Lichtheimia ramosa); Absidia is a recognised genus, so check the culture.', 'GmG_runs/reads/titration/n1500000/results.tsv; GmG_runs/external/zygolife/results.tsv'),
    ('A', 'lab_read_set+lab_assembly', 'NRRL A-20010', 'Gongronella sp. NRRL A-20010', 'GCF_025201355.1', 'Halteromyces radiatus', '0.9963 reads; 0.9997 assembly', 'ok (both)', 'READS+ASM',
     'reads and assembly both have Halteromyces radiatus as best hit', 'Halteromyces radiatus', 'Relabel in the lab sheet.', 'GmG_runs/reads/titration/n1500000/results.tsv; GmG_runs/external/zygolife/results.tsv'),
    ('A', 'lab_read_set+lab_assembly', 'NRRL 2543', 'Rhizomucor pusillus NRRL 2543', 'GCA_022702545.1', 'Phycomyces nitens', '0.9898 reads; 0.9986 assembly', 'ok (both)', 'READS+ASM',
     'reads and assembly both have Phycomyces nitens as best hit (different family from Rhizomucor)', 'Phycomyces nitens', 'Check the strain and relabel in the lab sheet.', 'GmG_runs/reads/titration/n1500000/results.tsv; GmG_runs/external/zygolife/results.tsv'),
    ('A', 'lab_read_set+lab_assembly', 'RSA 1136', 'Chaetocladium brefeldii RSA 1136', 'GCA_017871115.1', 'Cunninghamella blakesleeana', '0.9864 reads; 0.9892 assembly', 'ok (both)', 'READS+ASM',
     'reads and assembly both have Cunninghamella blakesleeana as best hit (different family)', 'Cunninghamella blakesleeana', 'Check the strain and relabel in the lab sheet.', 'GmG_runs/reads/titration/n1500000/results.tsv; GmG_runs/external/zygolife/results.tsv'),
    ('A', 'lab_read_set', 'JEL0774', 'Nowakowskiella sp. JEL0774', 'GCA_027604505.1', 'Gonapodya sp. JEL0774 (Monoblepharomycota)', '0.9903', 'ok; 3214 shared hashes', 'READS+STRAIN',
     'reads labelled Chytridiomycota (Nowakowskiella) match, at ANI 0.99, a reference genome that carries the SAME strain designation JEL0774 but the name Gonapodya (Monoblepharomycota); no ZyGoLife assembly of this strain was in the validation set', 'Gonapodya (strain JEL0774)', 'Check which name the culture JEL0774 carries; make the lab sheet and the NCBI genome agree.', 'GmG_runs/reads/titration/n1500000/results.tsv; reports/validation/reads_titration.md'),
    ('B', 'lab_read_set', 'JEL0847', 'Paraphysoderma sedebokerense JEL0847', 'GCA_027604345.1', 'Boothiomyces sp. JEL0866 (Chytridiomycota)', '0.9840', 'ok; 1332 shared hashes', 'READS',
     'reads only (no assembly in the set); best hit is another JEL strain, JEL0866, of a different genus', 'Boothiomyces', 'Check the culture identity (possible strain mix-up between JEL0847 and JEL0866).', 'GmG_runs/reads/titration/n1500000/results.tsv; reports/validation/reads_titration.md'),
    ('B', 'lab_read_set+lab_assembly', 'NRRL 3022', 'Conidiobolus incongruus NRRL 3022', 'GCA_054601915.1', 'Gilbertella persicaria (Mucoromycota)', '0.9332 reads; 0.9391 assembly', 'ok (both)', 'READS+ASM',
     'reads and assembly both put an Entomophthoromycota-labelled strain next to a Mucorales genome (different phylum) at ANI 0.93-0.94; ANI is below 0.98, so mislabel and Mucorales contamination are both possible', 'Gilbertella (Mucorales) or contaminated', 'Check culture purity and identity.', 'GmG_runs/reads/titration/n1500000/results.tsv; GmG_runs/external/zygolife/results.tsv'),
    ('C', 'lab_read_set', 'NRRL A-11121', 'Mucor mucedo NRRL A-11121', 'GCA_001986395.2', 'Irpex lacteus (Basidiomycota)', '0.9390', 'ok; 161 shared of 604 hashes (small sketch)', 'READS',
     'reads only; a Mucorales label gets a Basidiomycota best hit; few hashes, so likely a contaminated or mixed read set rather than a relabel', 'unknown (Irpex-like contaminant)', 'Check the read set for contamination before any relabel.', 'GmG_runs/reads/titration/n1500000/results.tsv; reports/validation/reads_titration.md'),
    ('C', 'lab_read_set', 'CBS 441.59', 'Syncephalastrum racemosum CBS 441.59', 'GCA_017871115.1', 'Cunninghamella blakesleeana (reads); Syncephalastrum racemosum (assembly)', '0.9736 reads; 0.9654 assembly', 'ok (both)', 'READS vs ASM',
     'the reads and the assembly of the same strain DISAGREE (the one case of 8 where they do): the assembly agrees with the label, the reads point at Cunninghamella', '', 'Check the read set (sample swap or contamination); the assembly is consistent with the label.', 'GmG_runs/reads/titration/n1500000/results.tsv; GmG_runs/external/zygolife/results.tsv'),
    ('C', 'lab_read_set+lab_assembly', 'NRRL 22899', 'Circinomucor circinelloides NRRL 22899 (reads sheet) / Mucor circinelloides (assembly sheet)', 'GCA_000401635.1', 'Mucor circinelloides', '0.9969 reads; 0.9992 assembly', 'ok (both)', 'READS+ASM',
     'the two lab sheets give two names for one strain; both calls are Mucor circinelloides (taxonomy: Circinomucor circinelloides is a newer name)', 'Mucor circinelloides (or Circinomucor circinelloides)', 'Keep; make the two lab sheets use one name.', 'GmG_runs/reads/titration/n1500000/results.tsv; GmG_runs/external/zygolife/results.tsv'),
    ('B', 'lab_assembly', 'NRRL 13127', 'Umbelopsis ovata NRRL 13127', 'GCA_054127015.1', 'Rhizopus delemar (Mucoromycetes; label is Umbelopsidomycetes)', '0.9976', 'ok', 'ASM',
     'one assembly; class differs from the label at ANI 0.9976 (a different class cannot be this similar)', 'Rhizopus delemar', 'Check the strain and relabel or flag as contaminated.', 'GmG_runs/external/zygolife/results.tsv; reports/validation/external_zygolife.md'),
    ('B', 'lab_assembly', 'NRRL 5842', 'Linnemannia nantahalensis NRRL 5842', 'GCA_027405735.1', 'Mucor racemosus (Mucoromycota; label is Mortierellomycota)', '0.9631', 'ok', 'ASM',
     'one assembly; phylum differs from the label at ANI 0.963', 'Mucor', 'Check the strain and relabel or flag as contaminated.', 'GmG_runs/external/zygolife/results.tsv'),
    ('B', 'lab_assembly', 'STR0027683', 'Trichoderma sp. TK1246 (IndigoAg, provider 2017 name)', 'GCA_060370295.1', 'Aspergillus fumigatus', '0.9893', 'ok', 'ASM',
     'one assembly; class differs (Sordariomycetes vs Eurotiomycetes) at ANI 0.989', 'Aspergillus fumigatus', 'Check the strain (provider name from 2017) and relabel.', 'GmG_runs/external/indigo/results.tsv, score.tsv'),
    ('B', 'lab_assembly', 'STR0027726', 'Alternaria sp. DC477 (IndigoAg)', 'GCA_023625515.1', 'Aspergillus terreus', '0.9864', 'ok', 'ASM',
     'one assembly; class differs (Dothideomycetes vs Eurotiomycetes) at ANI 0.986', 'Aspergillus terreus', 'Check the strain and relabel.', 'GmG_runs/external/indigo/results.tsv, score.tsv'),
    ('B', 'lab_assembly', 'STR0027671', 'Lophiostoma (IndigoAg)', 'GCA_028829675.1', 'Penicillium tannophilum', '0.9864', 'ok', 'ASM',
     'one assembly; class differs (Dothideomycetes vs Eurotiomycetes) at ANI 0.986', 'Penicillium tannophilum', 'Check the strain and relabel.', 'GmG_runs/external/indigo/results.tsv, score.tsv'),
    ('B', 'lab_assembly', 'STR0027783', 'Magnaporthe sp. QU1123 (IndigoAg)', 'GCA_001599375.1', 'Phoma herbarum', '0.9921', 'ok', 'ASM',
     'one assembly; class differs (Sordariomycetes vs Dothideomycetes) at ANI 0.992', 'Phoma herbarum', 'Check the strain and relabel.', 'GmG_runs/external/indigo/results.tsv, score.tsv'),
    ('B', 'lab_assembly', 'STR0027657', 'Magnaporthe sp. P3064 (IndigoAg)', 'GCA_027942835.1', 'Colletotrichum tropicale', '0.9884', 'ok', 'ASM',
     'one assembly; order differs (Magnaporthales vs Glomerellales, same class) at ANI 0.988', 'Colletotrichum tropicale', 'Check the strain and relabel.', 'GmG_runs/external/indigo/results.tsv, score.tsv'),
    ('C', 'lab_assembly', 'STR0027763', 'Neurospora sp. PS0169 (IndigoAg)', 'GCA_023626415.1', 'Penicillium sp. S1701', '0.9714', 'ok', 'ASM',
     'one assembly; class differs at ANI 0.971 (below 0.98)', 'Penicillium', 'Check if resources allow.', 'GmG_runs/external/indigo/results.tsv, score.tsv'),
    ('C', 'lab_assembly', 'STR0027688', 'Neosetophoma samarorum (IndigoAg)', 'GCA_023628855.1', 'Coniochaeta sp. PK30-01', '0.9781', 'ok', 'ASM',
     'one assembly; class differs (Dothideomycetes vs Sordariomycetes) at ANI 0.978; the hit itself is an unnamed Coniochaeta sp.', 'Coniochaeta', 'Check if resources allow.', 'GmG_runs/external/indigo/results.tsv, score.tsv'),
    ('C', 'lab_assembly', 'STR0027680', 'Diaporthe sp. P3581 (IndigoAg)', 'GCA_022385635.1', 'Xylaria scruposa', '0.9623', 'ok', 'ASM',
     'one assembly; order differs (same class) at ANI 0.962', 'Xylaria', 'Check if resources allow.', 'GmG_runs/external/indigo/results.tsv, score.tsv'),
    ('C', 'lab_assembly', 'STR0027681', 'Nemania sp. QU1039 (IndigoAg)', 'GCA_001879275.1', 'Coniochaeta ligniaria', '0.9548', 'ok', 'ASM',
     'one assembly; order differs (same class) at ANI 0.955', 'Coniochaeta', 'Check if resources allow.', 'GmG_runs/external/indigo/results.tsv, score.tsv'),
    ('C', 'lab_assembly', 'Cryphonectria__CPWV1912', 'Cryphonectria parasitica CPWV1912 (KassonLab)', 'n/a', 'Neonectria punicea', '0.9972', 'ok; input is archived herbarium DNA, a metagenome', 'ASM',
     'one assembly of herbarium metagenomic DNA: a Neonectria best hit is more likely a co-sequenced organism than a wrong name', '', 'Keep; treat as a mixed input.', 'GmG_runs/external/kasson/results.tsv; reports/validation/external_kasson.md'),
    ('C', 'lab_assembly', 'OphiocordycepsSymbiont__HPPF3', 'Ophiocordyceps sp. HPPF3 (KassonLab)', 'n/a', 'Massospora platypediae (Entomophthoromycota)', '0.9580', 'ok; fragmented; input is a metagenome', 'ASM',
     'one fragmented metagenome assembly; phylum differs from the label at ANI 0.958', '', 'Keep; treat as a mixed input.', 'GmG_runs/external/kasson/results.tsv; reports/validation/external_kasson.md'),
]
for t, kind, sample, lab, bacc, btax, ani, flags, et, evd, sug, act, src in LAB:
    add(tier=t, kind=kind, accession=f'(lab id) {sample}', bfd_id='(not in BFD)', current_organism=lab, best_hit_accession=bacc, best_hit_taxon=btax, ani=ani, flags=flags,
        evidence_type=et, evidence=evd, suggested_taxon=sug, suggested_action=act, **{'in_reference_2026.10.3': 'n/a (not in BFD)'}, sources=src)

# order: tier, kind, group
ordk = {'bfd_genome': 0, 'sra_run+bfd_genome': 1, 'lab_read_set+lab_assembly': 2, 'lab_read_set': 2, 'lab_assembly': 2, 'decoy_assembly': 3}
rows.sort(key=lambda r: (r['tier'], ordk[r['kind']], 0 if 'audit43' in r['flags'] else 1, -(float(re.match(r'[0-9.]+', r['ani']).group(0)) if r['ani'] and re.match(r'[0-9.]+', r['ani']) else 0)))
cnt = collections.Counter()
for r in rows:
    cnt[r['tier']] += 1
    r['row_id'] = f"{r['tier']}{cnt[r['tier']]:03d}"

with open(os.path.join(OUT, 'upstream_report.tsv'), 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=COLS, delimiter='\t', lineterminator='\n')
    w.writeheader()
    w.writerows(rows)

# ---------------------------------------------------------------- markdown tables
md = []
def esc(s):
    return str(s).replace('|', '/').replace('\n', ' ')

def table(sel, cols):
    out = ['| ' + ' | '.join(h for h, _ in cols) + ' |', '|' + '---|' * len(cols)]
    for r in sel:
        out.append('| ' + ' | '.join(esc(fn(r)) for _, fn in cols) + ' |')
    return out

short = lambda s, n=230: (s if len(s) <= n else s[:n - 3] + '...')
CM = [('id', lambda r: r['row_id']), ('accession', lambda r: r['accession']), ('current organism (BFD/NCBI)', lambda r: r['current_organism']),
      ('GmG best hit', lambda r: (r['best_hit_accession'] + ' ' + r['best_hit_taxon']).strip()), ('ANI', lambda r: r['ani']),
      ('size', lambda r: r['size_ratio']), ('flags', lambda r: r['flags']), ('evidence', lambda r: r['evidence_type']),
      ('suggested', lambda r: r['suggested_taxon'])]
def section(tier, title, kinds):
    sel = [r for r in rows if r['tier'] == tier and r['kind'] in kinds]
    md.append(f'#### {title} ({len(sel)})\n')
    if sel:
        md.extend(table(sel, CM))
        md.append('')
for t in 'ABC':
    md.append(f'### Tier {t}\n')
    section(t, 'BFD genomes (57 left out of the reference; the others kept, flag kept_in_reference or size_outlier)', {'bfd_genome'})
    section(t, 'Lab or provider labels of read sets and assemblies, the SRA run, the decoy', {'sra_run+bfd_genome', 'lab_read_set+lab_assembly', 'lab_read_set', 'lab_assembly', 'decoy_assembly'})
text = '\n'.join(md)


def grp(r):
    f = r['flags'].split(';')[0]
    if r['kind'] == 'bfd_genome':
        return {'audit43': 'BFD genomes left out: label audit (43)', 'truthset14': 'BFD genomes left out: truth-set run (14)',
                'kept_in_reference': 'BFD genomes kept in the reference: conflict with another kept genome', 'size_outlier': 'BFD genomes kept: size outliers only'}[f]
    if r['kind'].startswith('sra'):
        return 'SRA run (ERR14883901, assembly of the same BioSample)'
    if r['kind'] == 'decoy_assembly':
        return 'Decoy assembly GCF_965202585.1'
    return 'Lab or provider labels of read sets / assemblies (not NCBI accessions)'


GR = collections.OrderedDict()
for r in rows:
    GR.setdefault(grp(r), collections.Counter())[r['tier']] += 1
ctab = ['| group | A | B | C | total |', '|---|---|---|---|---|']
for k, c in GR.items():
    ctab.append(f"| {k} | {c['A']} | {c['B']} | {c['C']} | {sum(c.values())} |")
tot = collections.Counter(r['tier'] for r in rows)
ctab.append(f"| **all rows** | **{tot['A']}** | **{tot['B']}** | **{tot['C']}** | **{len(rows)}** |")
counts_text = '\n'.join(ctab)
p = os.path.join(OUT, 'upstream_report.md')
if os.path.exists(p):
    s = open(p).read()
    s = re.sub(r'<!-- TABLES -->.*<!-- /TABLES -->', lambda m_: '<!-- TABLES -->\n' + text + '\n<!-- /TABLES -->', s, flags=re.S)
    s = re.sub(r'<!-- COUNTS -->.*<!-- /COUNTS -->', lambda m_: '<!-- COUNTS -->\n' + counts_text + '\n<!-- /COUNTS -->', s, flags=re.S)
    open(p, 'w').write(s)

# counts
by = collections.Counter((r['tier'], r['kind'].split('+')[0] if r['kind'].startswith('bfd') else 'other', r['flags'].split(';')[0] if r['kind'] == 'bfd_genome' else '') for r in rows)
print(dict(cnt))
for k, v in sorted(by.items()):
    print(k, v)
print('total', len(rows))
