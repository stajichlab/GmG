#!/usr/bin/env python3
"""Pick ~100 Illumina WGS runs of BFD genomes with a trustworthy species label for the real-reads validation (#42).

usage: sra_plan.py ENA_RUNS.tsv ASSEMBLY_SUMMARY.txt INDEX_META.tsv EXCLUSIONS.txt OUT.tsv [N=100] [PER_GENUS=5] [PER_CENTER=25] [SEED=4]
  ENA_RUNS.tsv         ENA portal API `read_run` for tax_tree(4751), WGS, GENOMIC, ILLUMINA, PAIRED (fields in FIELDS below), tab separated
  ASSEMBLY_SUMMARY.txt NCBI genomes/genbank/fungi/assembly_summary.txt (BioSample of each assembly, organism name, genome size)
  INDEX_META.tsv       meta_full.tsv of the index release (id, lineage, species); EXCLUSIONS.txt the ids left out of the reference
A run is kept when its BioSample is the BioSample of a genome in the reference (matched by the assembly accession in the genome id), the species is
named the same way by three independent records (ENA scientific_name of the run, NCBI organism name of the assembly, the reference label), the name
is a binomial (not "sp."), the genome is not on the exclusion list, the run is paired, has >= 600,000 read pairs and a mean read length >= 90, and the
two FASTQ files are listed. One run per genome. Selection: genera in random order, at most PER_GENUS runs per genus, at most PER_CENTER per submitting
center, until N runs. Writes OUT.tsv (the accession list, committed in data/validation/sra_runs.tsv) and prints the spread (genera, phyla, centers).
"""
import csv, random, re, sys, collections

ena_p, asm_p, meta_p, excl_p, out = sys.argv[1:6]
N = int(sys.argv[6]) if len(sys.argv) > 6 else 100
PER_GENUS = int(sys.argv[7]) if len(sys.argv) > 7 else 5
PER_CENTER = int(sys.argv[8]) if len(sys.argv) > 8 else 25
rng = random.Random(int(sys.argv[9]) if len(sys.argv) > 9 else 4)
binom = lambda s: " ".join(re.sub(r"[\[\]]", "", s).split()[:2])


def good_name(s):
    w = binom(s).split()
    return len(w) == 2 and w[0][:1].isupper() and w[1][:1].islower() and w[1] not in ("sp.", "sp", "cf.", "aff.") and "." not in w[1]


excl = {l.strip() for l in open(excl_p) if l.strip()}
meta = {}
by_species = collections.Counter(); by_genus = collections.Counter()
for r in csv.DictReader(open(meta_p), delimiter="\t"):
    meta[r["id"]] = r
    by_species[r["species"]] += 1; by_genus[r["genus"]] += 1
# NCBI assembly summary: header is the second line ("#assembly_accession ...")
lines = open(asm_p).read().splitlines()
hi = next(i for i, l in enumerate(lines) if l.startswith("#assembly_accession"))
head = lines[hi].lstrip("#").split("\t")
asm = {}
for row in csv.reader(lines[hi + 1:], delimiter="\t"):
    d = dict(zip(head, row))
    if d.get("biosample"): asm[d["assembly_accession"]] = d
# reference genomes by BioSample
by_bs = collections.defaultdict(list)
for gid, m in meta.items():
    acc = re.match(r"^(GC[AF]_\d+\.\d+)", gid)
    if not acc or gid in excl or acc.group(1) not in asm: continue
    a = asm[acc.group(1)]
    if not (good_name(m["species"]) and good_name(a["organism_name"])): continue
    if binom(m["species"]) != binom(a["organism_name"]): continue
    by_bs[a["biosample"]].append((gid, a))
cands = {}
seen_genome = set()
for r in csv.DictReader(open(ena_p), delimiter="\t"):
    ids = {r["sample_accession"], r["secondary_sample_accession"]}
    hit = next((by_bs[b] for b in ids if b in by_bs), None)
    if not hit: continue
    gid, a = hit[0]
    m = meta[gid]
    if not good_name(r["scientific_name"]) or binom(r["scientific_name"]) != binom(m["species"]): continue
    try:
        rc, bc = int(r["read_count"]), int(r["base_count"])
    except ValueError:
        continue
    files = [f for f in r["fastq_ftp"].split(";") if f]
    f1 = next((f for f in files if f.endswith("_1.fastq.gz")), None); f2 = next((f for f in files if f.endswith("_2.fastq.gz")), None)
    if not (f1 and f2) or rc < 600_000 or bc / rc / 2 < 90: continue
    if gid in seen_genome: continue
    seen_genome.add(gid)
    gs = int(a.get("genome_size") or 0)
    cands[r["run_accession"]] = dict(run=r["run_accession"], biosample=a["biosample"], bfd_id=gid, phylum=m["phylum"], cls=m["class"], order=m["order"],
                                    family=m["family"], genus=m["genus"], species=m["species"], center=r["center_name"].strip() or "?",
                                    model=r["instrument_model"], read_pairs=rc, read_len=round(bc / rc / 2), genome_size=gs,
                                    species_others_in_ref=by_species[m["species"]] - 1, genus_others_in_ref=by_genus[m["genus"]] - 1, fastq_1=f1, fastq_2=f2)
print(len(cands), "eligible runs;", len({c["genus"] for c in cands.values()}), "genera,", len({c["phylum"] for c in cands.values()}), "phyla,",
      len({c["center"] for c in cands.values()}), "centers", file=sys.stderr)
genus_runs = collections.defaultdict(list)
for c in cands.values(): genus_runs[c["genus"]].append(c)
genera = sorted(genus_runs); rng.shuffle(genera)
pick, per_center = [], collections.Counter()
for rnd in range(PER_GENUS):
    for g in genera:
        lst = genus_runs[g]
        rng.shuffle(lst)
        if rnd < len(lst) and len(pick) < N and per_center[lst[rnd]["center"]] < PER_CENTER:
            pick.append(lst[rnd]); per_center[lst[rnd]["center"]] += 1
cols = ["run", "biosample", "bfd_id", "phylum", "cls", "order", "family", "genus", "species", "center", "model", "read_pairs", "read_len", "genome_size",
        "species_others_in_ref", "genus_others_in_ref", "fastq_1", "fastq_2"]
with open(out, "w") as fh:
    w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(pick)
print(len(pick), "runs;", len({c["genus"] for c in pick}), "genera;", dict(collections.Counter(c["phylum"] for c in pick)), ";",
      len({c["center"] for c in pick}), "centers:", dict(per_center.most_common(8)), file=sys.stderr)
