#!/usr/bin/env python3
"""Resolve the golden panel (#49) to NCBI FTP URLs and pin each file with its SHA-256.

usage: make_golden_panel.py OUT.tsv   (reads the ACCESSIONS table below; needs internet)
GenBank/RefSeq assembly files do not change for a given accession version, so URL + checksum pin the input of the regression test.
"""
import hashlib, re, sys, urllib.request

# accession, name, why it is in the panel
PANEL = [
    ("GCA_059919235.1", "Saccharomyces cerevisiae (IU_Scer_DoF1)", "species call, budding yeast, not in the reference"),
    ("GCA_056321735.1", "Schizosaccharomyces pombe", "species call, a different Ascomycota class"),
    ("GCA_051943995.1", "Malassezia restricta", "species call, Basidiomycota"),
    ("GCA_982436685.1", "Wallemia ichthyophaga", "species call, small Basidiomycota lineage"),
    ("GCA_985809295.1", "Vavraia culicis", "species call, Microsporidia"),
    ("GCF_000005845.2", "Escherichia coli K-12", "bacterium: no_match, no fungal call"),
    ("GCA_062046625.1", "Candidozyma auris", "species call, pathogen"),
    ("GCA_980652395.1", "Yamadazyma diddensiae", "species call, Pichiomycetes"),
    ("GCA_059056815.1", "Trichosporon austroamericanum", "species call, Tremellomycetes"),
    ("GCA_056632795.1", "Cryptococcus neoformans", "species call, pathogen"),
    ("GCA_060419215.1", "Rhodotorula sphaerocarpa", "species call, Microbotryomycetes"),
    ("GCA_977110955.1", "Rhizopus microsporus", "species call, Mucoromycota (thin lineage)"),
    ("GCA_985568295.1", "Vairimorpha sp.", "Microsporidia, provisional label"),
    ("GCA_062024532.1", "Entomophthoromycota assembly, distant from every genome", "distant decoy hit: the call is made without decoys, flag distant_non_fungal_hit"),
    ("GCA_061476465.1", "Ascomycota assembly with a low-confidence order call", "low_confidence: higher ranks only"),
    ("GCF_000146045.2", "Saccharomyces cerevisiae S288C (R64)", "in the reference: identical_to_reference"),
    ("GCF_000002765.5", "Plasmodium falciparum 3D7", "a decoy: not_fungal, identical_to_reference"),
]


def ftp_url(acc):
    p = f"{acc[:3]}/{acc[4:7]}/{acc[7:10]}/{acc[10:13]}"
    base = f"https://ftp.ncbi.nlm.nih.gov/genomes/all/{p}/"
    page = urllib.request.urlopen(base, timeout=60).read().decode()
    d = re.search(rf'href="({acc}_[^"/]+)/"', page).group(1)
    return f"{base}{d}/{d}_genomic.fna.gz"


out = sys.argv[1]
with open(out, "w") as f:
    f.write("id\tname\twhy\turl\tsha256\tbytes\n")
    for acc, name, why in PANEL:
        url = ftp_url(acc)
        data = urllib.request.urlopen(url, timeout=600).read()
        f.write(f"{acc}\t{name}\t{why}\t{url}\t{hashlib.sha256(data).hexdigest()}\t{len(data)}\n")
        print(acc, len(data), file=sys.stderr)
