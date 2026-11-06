// ENA metadata lookups for the popsketch page: turn run accessions into a label table (country, study, host, strain, species) from the
// ENA portal API, which allows cross-origin requests. Only metadata is fetched; the reads themselves are always the user's own files
// (streaming reads from ENA was tried and parked: its file server was too slow and unreliable to use from a browser).
// No network at import time; `fetchJson` is injectable so the functions can be tested without it.

const API = "https://www.ebi.ac.uk/ena/portal/api/search";
// below this a run is nearly empty (a few test reads), whatever else is true of it
const MIN_USABLE_READS = 100_000;
export const FIELDS = [
  "run_accession", "sample_accession", "study_accession", "scientific_name", "tax_id", "instrument_platform", "library_layout",
  "library_strategy", "library_source", "read_count", "base_count", "fastq_ftp", "fastq_bytes", "country", "strain", "isolate",
  "host", "collection_date", "center_name", "study_title",
];

// One token of the search box -> an ENA query clause. Run, study, sample accessions, a taxon id, or a species name.
export function clause(tok) {
  const t = tok.trim();
  if (/^[SED]RR\d+$/i.test(t)) return { q: `run_accession="${t.toUpperCase()}"`, taxon: false };
  if (/^(PRJ[EDN][A-Z]\d+|[SED]RP\d+)$/i.test(t)) return { q: `study_accession="${t.toUpperCase()}"`, taxon: false };
  if (/^(SAM[EDN][A-Z]?\d+|[SED]RS\d+)$/i.test(t)) return { q: `sample_accession="${t.toUpperCase()}"`, taxon: false };
  if (/^\d+$/.test(t)) return { q: `tax_tree(${t})`, taxon: true };
  if (/^[A-Za-z][A-Za-z .'-]+$/.test(t)) return { q: `tax_name("${t.replace(/"/g, "")}")`, taxon: true };
  throw new Error(`cannot read "${t}": give run (SRR…), study (PRJNA…), sample (SAMN…) accessions, a taxon id, or a species name`);
}

// Query text for the box. Taxon searches are limited to paired Illumina whole-genome reads, because that is what popsketch can use;
// accessions are taken as given (a run that cannot be used is shown with the reason instead of being hidden).
export function buildQuery(text) {
  // commas, semicolons and new lines separate items; spaces separate accessions but stay inside a species name
  const isAcc = (w) => /^([SED]RR\d+|PRJ[EDN][A-Z]\d+|[SED]RP\d+|SAM[EDN][A-Z]?\d+|[SED]RS\d+|\d+)$/i.test(w);
  const toks = text.split(/[,;\n]+/).flatMap((p) => {
    const words = p.trim().split(/\s+/).filter(Boolean);
    return words.length && words.every(isAcc) ? words : words.length ? [words.join(" ")] : [];
  });
  if (!toks.length) throw new Error("nothing to search for");
  const cl = toks.map(clause);
  const q = cl.map((c) => c.q).join(" OR ");
  const taxon = cl.some((c) => c.taxon);
  return taxon
    ? `(${q}) AND library_strategy="WGS" AND library_source="GENOMIC" AND instrument_platform="ILLUMINA" AND library_layout="PAIRED"`
    : `(${q})`;
}

export async function fetchJson(url) {
  const r = await fetch(url);
  if (!r.ok) throw new Error(`ENA answered ${r.status}`);
  const t = await r.text();
  return t.trim() ? JSON.parse(t) : [];
}

export async function findRuns(text, { limit = 2000, fetchJson: get = fetchJson } = {}) {
  const url = `${API}?result=read_run&query=${encodeURIComponent(buildQuery(text))}&fields=${FIELDS.join(",")}&limit=${limit}&format=json`;
  return annotate(await get(url));
}

// Parsed rows with numbers, the two FASTQ URLs, and why a run cannot be used (empty when it can).
export function annotate(rows) {
  return rows.map((r) => {
    const urls = (r.fastq_ftp || "").split(";").filter(Boolean).map((u) => (/^https?:/.test(u) ? u : `https://${u}`));
    const bytes = (r.fastq_bytes || "").split(";").filter(Boolean).map(Number);
    const paired = urls.filter((u) => /_[12]\.fastq\.gz$/.test(u));
    let why = "";
    if ((Number(r.read_count) || 0) < MIN_USABLE_READS) why = "fewer than 100,000 read pairs";
    else if (r.instrument_platform !== "ILLUMINA") why = `${r.instrument_platform || "unknown"} platform`;
    else if (r.library_strategy !== "WGS" || r.library_source !== "GENOMIC") why = `${r.library_strategy}/${r.library_source}, not genomic WGS`;
    else if (paired.length < 2) why = urls.length ? "single or unpaired FASTQ" : "no FASTQ at ENA";
    const mates = paired.length >= 2 ? [urls.find((u) => /_1\.fastq\.gz$/.test(u)), urls.find((u) => /_2\.fastq\.gz$/.test(u))] : [];
    return {
      ...r,
      read_count: Number(r.read_count) || 0,
      base_count: Number(r.base_count) || 0,
      urls: mates,
      mate_bytes: mates.map((m) => bytes[urls.indexOf(m)] || 0),
      why,
    };
  });
}

// Labels: id, country, study, host, strain, species (blank when ENA has nothing). Tabs and newlines in values would break the table.
const clean = (v) => String(v ?? "").replace(/[\t\r\n]+/g, " ").trim();
export function labelsText(runs) {
  const head = ["id", "country", "study", "host", "strain", "species"];
  const rows = runs.map((r) => [
    r.run_accession, clean(r.country).split(":")[0].trim(), clean(r.study_accession), clean(r.host), clean(r.strain || r.isolate), clean(r.scientific_name),
  ]);
  // a column that is blank for every sample would only be a label set with no groups
  const keep = head.map((_, i) => i === 0 || rows.some((r) => r[i]));
  const pick = (r) => r.filter((_, i) => keep[i]).join("\t");
  return [pick(head), ...rows.map(pick)].join("\n") + "\n";
}
