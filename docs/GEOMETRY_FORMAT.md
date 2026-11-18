# Fixed production geometry format

As of 2026-10-21, the default builder, reference-store API/CLI and production
examples use **packed, zlib level 9**, with the existing 65,536-unit chunks.
This follows the previous auto objective: smallest validated encoded payload.
Production builds evaluate only that one candidate. Full array equality,
read-back/replay verification, source identities and payload checksums remain.
No coverage, filtering, sampling-unit or numerical semantics change.

Evidence: the completed IH batch contains 26,065 sample–GCF records, 4,557
samples and 906 GCFs. All 26,065 auto winners were packed: level9 25,677;
level6 383; level1 5. Auto's stable tie order accounts for some non-level9 winners.
Fixed packed9 totals 37,193,285,516 bytes, only 134 bytes more than choosing the
minimum per sample. Candidate bytes include encoded numeric files and storage.json,
exclude common bundle metadata and do not estimate the entire library footprint.

| Candidate | Total encoded bytes | Recorded write seconds, summed |
|---|---:|---:|
| packed1 | 41,779,247,293 | 3,794.97 |
| packed6 | 37,867,112,376 | 5,160.58 |
| packed9 | 37,193,285,516 | 12,247.92 |
| NPZ6 | 102,533,099,248 | 78,253.98 |
| NPY0 | 542,561,721,835 | 2,268.18 |

Packed6 writes faster in these observations, with about 1.81% more bytes than
packed9. The chosen production policy prioritizes compact storage and removes
repeated format searches; it is not a claim that level9 is the fastest encoder.
Recorded candidate write/load times exclude replay and other overhead, overlap
across concurrent workers, and cannot predict end-to-end speedup.

`--geometry-format packed` (builder) and `--format packed` (reference_store)
now mean fixed packed9; both are defaults. Explicit `auto`, `npy` and `npz`
remain available for controlled diagnostics/compatibility, without silent fallback.
Auto still tests its original five candidates and picks the minimum by bytes.
Readers continue using storage.json and can read all previously saved formats
and compression levels. Existing libraries, frozen plans and outputs are not
rewritten. New source fingerprints require new construction plans; retain old
source for active workers. Nothing is deployed to an external /groups source
copy implicitly.
