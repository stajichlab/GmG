// Read the first `maxReads` records of a FASTQ file (plain or gzip) for popsketch: exactly that many, then stop reading the file.
// Used by popsketch-worker.js; a module of its own so that node can test it.
import { openDecoded } from "./gunzip-browser.js";

export const MAX_HEAD_BYTES = 1.5 * 2 ** 30; // bytes of reads held for one file: more is an error, not a silent truncation
export const MAX_LINE_BYTES = 64 * 2 ** 20; // a FASTQ line this long without a line break is not a FASTQ file

const NL = 10;

/**
 * The first `maxReads` FASTQ records of `file` (all of it when it is not FASTQ: a FASTA genome is read whole), as one Uint8Array.
 * - Exactly `maxReads` records: the last chunk is cut at the end of the record, so the count does not depend on the chunk size.
 * - Gzip is recognised by its first bytes. The CRC of a member is checked when the member is read to its end; a file cut at the read
 *   limit is not read to its end, so what follows the cut is not checked and an error that lies beyond it is never seen.
 * - Past `maxBytes` of reads the file is cut at the end of a record with a note (`onNote`); a non-FASTQ file that large is an error.
 *   A line of more than MAX_LINE_BYTES with no line break is an error.
 * `onNote(message)` also receives the decoder's warning about bytes after the last gzip member.
 */
export async function readHead(file, maxReads, { maxBytes = MAX_HEAD_BYTES, maxLine = MAX_LINE_BYTES, onNote = () => {} } = {}) {
  if (!(Number.isInteger(maxReads) && maxReads >= 1)) throw new Error(`reads per file must be a whole number of 1 or more, not ${maxReads}`);
  const stream = await openDecoded(file, { onWarn: onNote });
  const reader = stream.getReader();
  const maxLines = maxReads * 4;
  const chunks = [];
  let total = 0, lines = 0, lineLen = 0, recordEnd = 0, fastq = null;
  for (;;) {
    const { done, value } = await reader.read();
    if (done) break;
    if (fastq === null && value.length) fastq = value[0] === 0x40; // '@'
    if (fastq) {
      let start = 0;
      for (let i = 0; i < value.length; i++) {
        if (value[i] !== NL) continue;
        lines++;
        lineLen = 0;
        start = i + 1;
        if (lines % 4 === 0) recordEnd = total + i + 1; // offset just after the last complete record
        if (lines === maxLines) { // the end of the last wanted record: keep up to here and stop reading
          chunks.push(value.subarray(0, i + 1));
          await reader.cancel().catch(() => {}); // a stream that already failed (an error beyond the cut) is not our concern
          return join(chunks, total + i + 1);
        }
      }
      lineLen += value.length - start;
      if (lineLen > maxLine) throw new Error(`${file.name}: a line of more than ${maxLine / 2 ** 20} MiB with no line break: not a FASTQ file`);
    }
    chunks.push(value);
    total += value.length;
    if (total > maxBytes) {
      await reader.cancel().catch(() => {}); // a stream that already failed (an error beyond the cut) is not our concern
      if (!fastq) throw new Error(`${file.name}: more than ${(maxBytes / 2 ** 30).toFixed(1)} GiB of sequence: too large for the browser`);
      onNote(`${file.name}: stopped at ${Math.floor(lines / 4)} reads: the ${(maxBytes / 2 ** 30).toFixed(1)} GiB limit for one file was reached before ${maxReads} reads`);
      return join(chunks, total).subarray(0, recordEnd);
    }
  }
  if (total === 0) throw new Error(`${file.name}: empty file`);
  return join(chunks, total);
}

function join(chunks, length) {
  const bytes = new Uint8Array(length);
  let off = 0;
  for (const c of chunks) {
    if (off + c.length > length) { bytes.set(c.subarray(0, length - off), off); break; }
    bytes.set(c, off);
    off += c.length;
  }
  return bytes;
}

/** All of `file`, decoded when it is gzip, with the CRC of every gzip member checked (the whole file is read, so it can be). */
export async function readWhole(file, { onNote = () => {} } = {}) {
  const reader = (await openDecoded(file, { crc: true, onWarn: onNote })).getReader();
  const chunks = [];
  let total = 0;
  for (;;) {
    const { done, value } = await reader.read();
    if (done) break;
    chunks.push(value);
    total += value.length;
  }
  return join(chunks, total);
}
