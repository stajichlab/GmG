// The multi-member gunzip (gunzip.js) with pako (vendor/pako_inflate.mjs, the unmodified pako 1.0.11 inflate build wrapped as
// an ES module so that it behaves the same in a page, a worker and node).
import pako from "./vendor/pako_inflate.mjs";
import { makeGunzip } from "./gunzip.js";

export const gunzip = makeGunzip(pako);
