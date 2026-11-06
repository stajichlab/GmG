import { test } from "node:test";
import assert from "node:assert/strict";
import { buildQuery, annotate, labelsText, findRuns, clause } from "../../web/ena.js";

const row = (o) => ({
  run_accession: "SRR1", sample_accession: "SAMN1", study_accession: "PRJNA1", scientific_name: "Cercospora beticola", instrument_platform: "ILLUMINA",
  library_layout: "PAIRED", library_strategy: "WGS", library_source: "GENOMIC", read_count: "5000000", base_count: "1000000000",
  fastq_ftp: "ftp.sra.ebi.ac.uk/a/SRR1_1.fastq.gz;ftp.sra.ebi.ac.uk/a/SRR1_2.fastq.gz", fastq_bytes: "500000000;490000000",
  country: "USA:North Dakota", strain: "09-40", isolate: "", host: "Beta vulgaris", ...o,
});

test("accessions of every kind, a taxon id and a species name become one query", () => {
  assert.equal(buildQuery("SRR1713175"), '(run_accession="SRR1713175")');
  assert.equal(buildQuery("SRR1 SRR2 SRR3"), '(run_accession="SRR1" OR run_accession="SRR2" OR run_accession="SRR3")');
  assert.equal(buildQuery("PRJNA270309, SAMN03265455"), '(study_accession="PRJNA270309" OR sample_accession="SAMN03265455")');
  const q = buildQuery("5518");
  assert.match(q, /tax_tree\(5518\)/);
  assert.match(q, /library_layout="PAIRED"/);
  assert.match(buildQuery("Cercospora beticola"), /tax_name\("Cercospora beticola"\)/);
  assert.throws(() => buildQuery("   "), /nothing/);
  assert.throws(() => clause("SRR1; DROP"), /cannot read/);
});

test("runs get a reason when they cannot be used for reads", () => {
  const [ok, nano, tiny, other] = annotate([
    row({}),
    row({ instrument_platform: "OXFORD_NANOPORE" }),
    row({ read_count: "300" }),
    row({ library_strategy: "RNA-Seq", library_source: "TRANSCRIPTOMIC" }),
  ]);
  assert.equal(ok.why, "");
  assert.match(nano.why, /OXFORD_NANOPORE/);
  assert.match(tiny.why, /fewer than/);
  assert.match(other.why, /not genomic WGS/);
});

test("label table: country trimmed to its first part, empty columns dropped, no tabs inside values", () => {
  const rs = annotate([row({ run_accession: "A", strain: "x\ty" }), row({ run_accession: "B", country: "", host: "" })]);
  const lines = labelsText(rs).trim().split("\n");
  assert.equal(lines[0], "id\tcountry\tstudy\thost\tstrain\tspecies");
  assert.equal(lines[1], "A\tUSA\tPRJNA1\tBeta vulgaris\tx y\tCercospora beticola");
  assert.equal(lines[2], "B\t\tPRJNA1\t\t09-40\tCercospora beticola");
  const bare = labelsText(annotate([row({ country: "", host: "", strain: "", isolate: "" })])).trim().split("\n");
  assert.equal(bare[0], "id\tstudy\tspecies");
});

test("findRuns asks the ENA portal for the run table and parses what comes back", async () => {
  let seen = "";
  const runs = await findRuns("SRR1", { fetchJson: async (u) => { seen = u; return [row({})]; } });
  assert.match(seen, /^https:\/\/www\.ebi\.ac\.uk\/ena\/portal\/api\/search\?result=read_run&query=/);
  assert.match(seen, /fields=run_accession,sample_accession/);
  assert.equal(runs.length, 1);
  assert.equal(runs[0].read_count, 5000000);
});
