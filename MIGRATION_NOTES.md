# Migration Notes: `pattyfam_compute` (Perl) → `signature_kmers` (parallel C++)

Status report and task list for the ongoing port of the PATRIC protein-family
(`pattyfam`) pipeline from Perl into efficient parallel C++. Written 2026-07 from a
source-and-script audit of both modules. Intended as a working reference for later
migration work — verify line numbers against the current tree before acting.

Modules:
- `signature_kmers` (**SIG**) — `dev_container/modules/signature_kmers`. The
  consolidation target: modern C++20, cmph perfect-hash + NuDB backends, TBB-parallel.
- `pattyfam_compute` (**PF**) — `dev_container/modules/pattyfam_compute`. Being reduced
  to a Perl orchestration layer; still carries an older C++ generation.

---

## 1. Executive summary

The migration runs in one direction: **PF → SIG**. It is at two different stages:

- **Pipeline (Perl → binaries): essentially complete.** Every live `pf-*.pl` step
  shells out to a SIG binary. PF's own C++ binaries are compiled but invoked by no
  live script.
- **Source tree (duplicated C++ infrastructure): partial and diverged.** PF keeps its
  own, older, binary-incompatible copies of shared headers and two whole kmer engines.
  Some functionality (weighting, OTU, DNA six-frame, FASTQ) exists only in PF; one
  tool (`propagate_names`) is a C++ port that is compiled but not yet wired in.

The overall shape: Perl is becoming glue / I/O / SEED-DB marshaling; C++ absorbs the
compute-bound inner loops (kmer build, annotation, all-to-all distance).

---

## 2. Pipeline wiring (what the live Perl calls)

| Pipeline step | Script (file:line) | Binary | Module |
|---|---|---|---|
| Build kmer DB | `pf-compute-signature-kmers.pl:263` | `kmers-build-signatures` | SIG |
| Annotate / call | `pf-compute-local-families.pl:143` | `kmers-annotate-seqs` | SIG |
| Distance (folder) | `pf-compute-kmer-distances.pl:376` | `kmers-matrix-distance-folder` | SIG |
| Distance (merge) | `pf-merge-stage-1.pl:191`, `pf-merge-stage-1-guts.pl:127` | `kmers-matrix-distance-files` | SIG |

- Legacy SEED `kmer_search` is retired — only a commented-out line remains
  (`pf-annotate-seqs.pl:151`); `pf-annotate-seqs.pl` itself is the superseded
  server-based path, replaced by `pf-compute-local-families.pl → kmers-annotate-seqs`.
- PF binaries `build_signature_kmers`, `nudb_call_kmers`, `propagate_names`, and the
  top-level `recall_proteins` are **built but not invoked** by any live script.

---

## 3. Evidence of translation in progress (in `pattyfam_compute/scripts`)

1. **`.pl~` backups are literal pre-port snapshots.** `pf-compute-signature-kmers.pl~`
   selected PF-local binaries by run-mode (`build_signature_kmers` / `recall_proteins`);
   the live `.pl:263` calls `kmers-build-signatures` with the new
   `--perfect-hash`/`--perfect-hash-data`, `--ignored-functions-file`,
   `--function-override-file` options.
2. **Perl as a shim around un-ported C++.** `pf-compute-signature-kmers.pl:194-229`
   de-duplicates per-ID functions in Perl, with the comment: *"The C++ code currently
   does not handle that case and as such may compute some statistics incorrectly."* It
   also filters deleted FIDs and rewrites FASTA in Perl before handing clean inputs to
   the C++ (changed from `symlink` to explicit rewrite).
3. **Commented-out flags for not-yet-ported features**, e.g. `--truncated-pegs` in
   `pf-compute-local-families.pl:145` (truncated-genes concept exists in Perl —
   `pf-find-truncated-genes.pl` — but is not yet wired into the C++ annotate call).
4. **Hot loops pushed into C++.** Merge-stage distance streams a file list into
   `kmers-matrix-distance-files` via `/dev/fd/1` piped to `mcl`
   (`pf-merge-stage-1-guts.pl:110-135`). Git milestones: PF `97f455d "Inner loop
   compute for merge stages 1 & 2"`; SIG "New method for computing distance without
   O(n^2) step".
5. **Whole-tool port.** `pattyfam_compute/src/propagate_names.cc` is the C++ of the
   Perl `p3x-propagate-family-names` (referenced by
   `pf-update-propagated-family-names.pl:2`) — compiled but not yet invoked.

---

## 4. Source-tree divergence (the hard part still to do)

### 4.1 Duplicated infrastructure headers
PF's `Makefile:39` adds `-I../signature_kmers/src`, but PF/src keeps its own copies of
several headers. Because same-directory `#include "x.h"` shadows the `-I` path, PF's
binaries compile against **PF's own copies**, and the include flag is effectively inert
for them.

| File | Status |
|---|---|
| `seed_utils.h` | identical |
| `operators.h` | lightly diverged (SIG made `split` inline) |
| `fasta_parser.h` / `.cc` | lightly/moderately diverged (SIG added a defline-parser state) |
| `function_map.h` | substantially diverged (SIG: per-function length stats, typed `FunctionIndex`=uint32 + overflow check, const API; PF: `unsigned short`) |
| `nudb_kmer_db.h` | substantially diverged (different template signature + stored record) |

### 4.2 Binary-incompatible stored record
| | SIG `StoredKmerData` (`kmer_data.h:114`) | PF `NuDBKmerDb::KData` (`pattyfam_compute/src/nudb_kmer_db.h`) |
|---|---|---|
| fields | avg_from_end, function_index(**u32**), mean, median, var | otu_index, avg_from_end, function_index(**u16**), **weight(float)**, mean, median, var, **n_proteins** |
| size | ~12–16 B | ~24 B |

Databases written by one tree cannot be read by the other. **Do not attempt to unify
by deleting PF's duplicate headers** — that switches resolution to SIG's incompatible
layout. Reconcile the record format and index widths first.

### 4.3 Two kmer-engine generations still in PF
- `kguts.{h,cc}` — classic SEED KmerGuts: base-20 encoded 8-mers, on-disk hash image,
  **DNA six-frame translation, OTU tracking, weighted scoring**.
- `kmer_generic.{h,tcc}` + `kmer_nudb.h` — transitional generic-over-Caller engine
  (what `nudb_call_kmers` uses).
- SIG collapsed all of this into one `FunctionCaller<KmerDb>` templated over
  `CmphKmerDb`/`NuDBKmerDb`, amino-acid FASTA only.

---

## 5. Feature deltas (SIG vs PF)

**Dropped by SIG (present only in PF):**
- Signature **weighting** (naive-Bayes log-odds; see `kmer_derivation.pdf` and §7).
- **OTU** tracking (`otu_index`, OTU stats).
- **DNA / six-frame** translation (KmerGuts only).
- **FASTQ** input (`fastq_parser`, compiled in PF but unused by its generic caller).

**Added by SIG (not in PF):**
- cmph **perfect-hash** backend + mmap flat data file.
- **Off-length** (short/long) call gating via protein-length median/MAD.
- **TBB-parallel** fasta processing + lock-free output queue.
- `kmers-matrix-distance*` suite (all-to-all distance, O(n²)-avoiding method).
- `--function-override-file`, per-function length statistics.

---

## 6. Open task list for finishing the migration

Roughly ordered; each is independent unless noted.

1. **Decide the fate of dropped features.** Confirm whether **weighting**, **OTU**, and
   **DNA/six-frame** are intentional removals or gaps the family pipeline still needs.
   - Weighting: good results reported without it (count-based scoring). If re-added,
     implement in SIG and **fix the occurrence-vs-sequence counting bug** — see
     `kmer_derivation.pdf` §8.2 and §7. Sketch: in `signature_build.tcc::process_kmer_set`
     count **distinct `seq_id`** for `NSi`/`NSiFj` (not occurrences); restore the two
     `KeptKmer` count fields (`signature_build.h:43-44`) as sequence counts; add
     `float weight` to `StoredKmerData` (breaks DB compat — version/rebuild); add a
     finalize `parallel_for` computing weight; add an opt-in weighted path in
     `call_functions`. The 80% prune's counting unit (§7.1) is a separate,
     behavior-changing decision — do the weight fix first, measure recall deltas.
2. **Port `propagate_names` into the live pipeline** (or into SIG). It is a compiled
   C++ port of `p3x-propagate-family-names` with no live caller; decide its home and
   wire `pf-update-propagated-family-names.pl` (or a successor) to it.
3. **Remove PF's un-ported C++ compensations from Perl.** The per-ID function dedup in
   `pf-compute-signature-kmers.pl:194-229` exists because the C++ mishandles repeated
   IDs. Fix `kmers-build-signatures` to handle duplicate assignments, then drop the
   Perl workaround.
4. **Wire deferred flags.** `--truncated-pegs` (annotate + distance) is stubbed out;
   implement in the C++ or delete the dead flag references.
5. **Reconcile / delete PF's duplicate infrastructure.** Once the stored record and
   index widths are unified, remove PF's shadowing copies of `function_map.h`,
   `nudb_kmer_db.h`, `fasta_parser.*`, `operators.h`, `seed_utils.h` so PF genuinely
   compiles against SIG headers (its `-I../signature_kmers/src` already intends this).
6. **Retire PF's dead binaries.** After 1–5, drop `build_signature_kmers`,
   `nudb_call_kmers`, and `recall_proteins` (and the `kguts`/`kmer_generic` engines)
   from PF's `Makefile:7` and delete the sources — unless a feature from §5 must be
   preserved from them first.
7. **Declare the dependency.** `pattyfam_compute/DEPENDENCIES` lists only `sys_tools`;
   add `signature_kmers` (build-time headers + runtime binaries on PATH).
8. **Rebuild databases.** Any change to `StoredKmerData` invalidates existing
   `kmer_data.{mph,dat}` / NuDB stores — plan a coordinated rebuild + version marker.

---

## 7. Key file references

SIG: `src/kmers-build-signatures.cc`, `signature_build.{h,tcc}`, `call_functions.{h,tcc}`,
`kmer_data.h` (`StoredKmerData` :114), `cmph_kmer.h`, `nudb_kmer_db.h`,
`matrix_distance.h`, `kmer_derivation.pdf`.

PF: `src/build_signature_nudb.cc` (active build; `compute_weight_of_signature` :764-776),
`src/build_signature_kmers.cc` (dead), `src/nudb_call_kmers_generic.cc` (active caller),
`recall_proteins.cc` (top-level, build+recall), `src/kguts.*`, `src/kmer_generic.*`,
`src/nudb_kmer_db.h` (`KData`), `src/propagate_names.cc`, `Makefile:7,39`,
`DEPENDENCIES`, `scripts/pf-compute-signature-kmers.pl`, `scripts/pf-compute-local-families.pl`,
`scripts/pf-compute-kmer-distances.pl`, `scripts/pf-merge-stage-1{,-guts}.pl`.

Background: `kmer_derivation.pdf` (weight derivation + known discrepancies);
see also this module's `CLAUDE.md` "Relationship to the original Perl build".
