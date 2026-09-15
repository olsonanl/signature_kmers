# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

`signature_kmers` is a C++20 module in the BV-BRC/PATRIC `dev_container`. It builds
and queries **signature kmers**: length-8 protein kmers that are diagnostic of a
specific protein function, used to annotate protein sequences and to compute
sequence-to-sequence distances. It is a performance-focused reimplementation of
older SEED/pattyfam kmer tooling (`kmers-annotate-seqs` is a drop-in for
`pf-annotate-seqs`).

## Build

This is a dev_container module, not standalone. The `Makefile` includes
`../../tools/Makefile.common`, which supplies `KB_RUNTIME`, `BOOST`, `TBB`, and
compiler variables. Before building you must source the environment from the
`dev_container` root:

```bash
source ../../user-env.sh      # sets KB_TOP, KB_RUNTIME, paths
make                          # builds all APP_CXX programs into ./ (and BIN_DIR)
make <program-name>           # build a single program, e.g. make kmers-call-functions
```

Build dependencies beyond Boost/TBB:
- **NuDB** — header-only; `make` clones it via `git clone` into `./NuDB` if absent.
- **cmph** (perfect hashing) — built locally from `cmph-2.0.2.tar.gz`; the resulting
  `lib/libcmph.*` and `include/cmph.h` are linked/included directly from this dir
  (`CMPH = $(shell pwd)` in the Makefile). If a link fails on missing cmph, rebuild
  it from the unpacked `cmph-2.0.2/` tree.

Compiler flags of note: `-std=c++20 -O3 -D_GLIBCXX_ASSERTIONS`, TBB deprecation
messages suppressed. There is no test target and no linter; `make clean` is a no-op.
Object/dependency files (`src/*.o`, `src/*.d`) are gitignored and `.d` files are
`include`d by the Makefile, so a clean checkout builds from scratch.

Ad-hoc scratch programs (`x.cc`, `taccum.cc`, `src/test-fasta-parse.cc`, etc.) exist
at the repo root and are not part of `APP_CXX` — ignore them unless asked.

## Programs

Perl (in `scripts/`, wrapped and deployed by `make`):
- `kmers-setup-build.pl` — runs **inside the SEED environment** (uses `FIG`,
  `FIG_Config`). Extracts per-genome `Seqs/` (fasta) and `Annotations/0/` files that
  become the training input for `kmers-build-signatures`. Not buildable/runnable
  outside a SEED install.

C++ (`src/*.cc`, one `main` each):
- `kmers-build-signatures` — the training step. Consumes function-definition files +
  fasta, produces the kmer database and metadata (see Data flow).
- `kmers-call-functions` — annotate fasta sequences with function calls from a built
  database.
- `kmers-annotate-seqs` — pattyfam-pipeline annotation entry point (per-genus dirs,
  offlength handling); drop-in for `pf-annotate-seqs`.
- `kmers-matrix-distance` / `-files` / `-folder` / `-merge` — all-to-all protein
  distance via shared signature-kmer counts. The current distance method (kmer →
  set-of-sequences inversion, avoiding the O(n²) pairwise scan) lives in
  `kmers-matrix-distance.cc` and the reusable `src/matrix_distance.h`.

## Architecture

The design pivots on two interchangeable kmer-database backends that expose the same
`fetch`/`insert` surface, so `FunctionCaller` and the build code are templated over
the DB type:

- **`NuDBKmerDb<StoredData,K>`** (`nudb_kmer_db.h`) — NuDB key-value store. Used as a
  **write** target during the build (`--nudb-file`).
- **`CmphKmerDb<StoredData,K>`** (`cmph_kmer.h`) — a cmph **perfect hash** (`.mph`)
  plus a flat, `mmap`-ed data file (`.dat`) of `StoredKmerData` indexed by hash slot.
  This is the **read/query** backend used by all the calling/distance tools. It
  `madvise(MADV_POPULATE_READ)`s the mapping for fast lookups. `kmers-call-functions`,
  `-annotate-seqs`, and `-matrix-distance` all open `data_dir/kmer_data.{mph,dat}`.

Core types (`kmer_data.h`):
- `Kmer<K> = std::array<char,K>`; `K` is fixed at **8** throughout the programs.
- `FunctionIndex` (uint32) indexes into `function.index`; `UndefinedFunction` is the
  max sentinel. `StoredKmerData` is the 12-byte per-kmer record actually persisted
  (function index, avg-from-end, mean/median/variance of protein length).
- `for_each_kmer<N>` walks a sequence's kmers, skipping windows containing `X`/`*`.

Build pipeline — `SignatureBuilder<K>` (`signature_build.h` + `signature_build.tcc`),
driven by `kmers-build-signatures.cc`:
1. `load_function_data` — build the id→function map and function→genome evidence via
   `FunctionMap` (`function_map.h`).
2. `load_fasta` — read training fasta (regular set + "keep-functions" set).
3. `process_kept_functions` — decide which functions have enough genome
   representation (`--min-reps-required`) to earn signatures; writes `function.index`.
4. `extract_kmers` / `process_kmers` — collect kmer→attributes, keep only kmers
   diagnostic of a single function (the signatures), compute per-kmer stats.
5. Outputs into `--kmer-data-dir`: `function.index`, `distinct_functions`,
   `recall.report.d/` (self-recall diff of calls vs. original annotations),
   `final.kmers`, optionally the NuDB db, and the cmph perfect hash + data file
   (`perfect_hash.h::build_perfect_hash`).

Calling pipeline — `FunctionCaller<KmerDb>` (`call_functions.h` +
`call_functions.tcc`): streams fasta (`fasta_parser.h/.cc`), for each sequence
gathers kmer hits (`hit_cb`), groups them into `KmerCall` runs, and `find_best_call`
picks the winning function with a score. Callers pass a `hit_cb` and a `call_cb`;
output is assembled off-thread and drained by a single writer thread through a
`tbb::concurrent_bounded_queue` of `boost::asio::streambuf` buffers.

Concurrency: the whole module leans hard on Intel **TBB** — `tbb::parallel_for`
over file/kmer ranges, `concurrent_unordered_map/set/vector`, and
`tbb::global_control` to cap threads (`-j`/`--n-threads`). Assume shared containers
are the concurrent variants and keep new code thread-safe accordingly.

Supporting headers: `seq_id_map.h` (int↔string id interning for the distance
matrix), `calc_natural_breaks.h` (Jenks breaks), `path_utils.h`
(`populate_path_list`, directory helpers), `operators.h`, `seed_utils.h`.

## Relationship to the original Perl build

This module is a port/reimplementation of the SEED Perl signature-kmer build. The
original toplevel is
`~olson/FIGdisk/dist/releases/current/FigKernelScripts/rebuild_converged_Data.pl`.
Understanding that pipeline explains several design choices here:

- The original build is a **convergence** process, not a single pass.
  `rebuild_converged_Data.pl` first materializes per-genome `Seqs/` (fasta) and
  `Annotations/0/` (assignments) from a chosen SEED, then runs
  `run1_convergence_step` three times to produce `Data.0 → Data.1 → Data.2`. Each
  step builds a kmer set, **recalls every genome against it** (`kmer_search`), merges
  the new calls into `Annotations/N`, and records changed calls in `New/N`. Feeding
  corrected annotations back in is the "convergence." Finally
  `compute_stats_for_kmer_hits` computes the length/density statistics used for
  z-scores.
- Port mapping:
  - `scripts/kmers-setup-build.pl` ports the **data-prep** portion of
    `rebuild_converged_Data.pl` (the `Seqs/` + `Annotations/0/` generation, fragment
    /frameshift suffixing, additional-fasta and function-override handling). It still
    requires a live SEED (`FIG`, `FIG_Config`).
  - `kmers-build-signatures` ports the kmer-construction step (`km_build_Data` /
    `kmer_search --allow-rebuild`) **plus** the self-recall. The `recall.report.d/`
    output and the `saver`/`recall_report` structs — which record IDs whose new call
    differs from the original stripped function — are the C++ equivalent of the
    per-step `New/N` "changed calls" files that drive convergence.
  - `kmers-call-functions` / `kmers-annotate-seqs` port `kmer_search` (annotation).
  - The mean/median/variance-of-protein-length fields in `StoredKmerData`, and the
    off-length hit filtering in the callers, are the C++ equivalent of
    `compute_stats_for_kmer_hits` (z-scores on protein length / hit density).

The C++ `kmers-build-signatures` is deliberately **single-pass**: it does one build
plus an in-process recall, and does not iterate to convergence. The convergence loop
existed because the kmer method performs best on *consistent* annotations, and the
recall step was what forced consistency across the training set. Since the project
began, the reference CoreSEED database has become much more consistent and accurate,
so the value of iterating is now unclear — treat the single-pass design as
intentional, not an unfinished port. The recall output (`recall.report.d/`) is still
produced and is useful for inspecting where calls disagree with input annotations.

Background reading:
- Edwards et al., "Real Time Metagenomics: Using k-mers to annotate metagenomes,"
  Bioinformatics (2012) — https://pmc.ncbi.nlm.nih.gov/articles/PMC3519453/
  (the k-mer annotation method itself).
- Overbeek et al., "The Subsystems Approach to Genome Annotation and its Use in the
  Project to Annotate 1000 Genomes," Nucleic Acids Res. 33(17):5691–5702 (2005) —
  https://pmc.ncbi.nlm.nih.gov/articles/PMC1251668/ (the subsystems/annotation
  framework the functions and roles come from).
- Overbeek et al., "The SEED and the Rapid Annotation of microbial genomes using
  Subsystems Technology (RAST)," Nucleic Acids Res. 42(Database issue):D206–D214
  (2014) — https://pmc.ncbi.nlm.nih.gov/articles/PMC3965101/ (SEED/RAST context).

## Related module: `../pattyfam_compute` (active migration)

`signature_kmers` is the consolidation target of an in-progress port of the PATRIC
protein-family (`pattyfam`) pipeline from Perl into efficient parallel C++.
`pattyfam_compute` is being reduced to a Perl orchestration layer that shells out to
the `kmers-*` binaries built here. Practical implications when working in either tree:

- The **live** pattyfam pipeline already calls this module's binaries:
  `pf-compute-signature-kmers.pl → kmers-build-signatures`,
  `pf-compute-local-families.pl → kmers-annotate-seqs`,
  `pf-compute-kmer-distances.pl → kmers-matrix-distance-folder`,
  `pf-merge-stage-1{,-guts}.pl → kmers-matrix-distance-files`.
- `pattyfam_compute/src` still contains an **older, diverged** generation of this code
  (`kguts.*`, `kmer_generic.*`, `build_signature_nudb.cc`, `nudb_call_kmers*`, and its
  own copies of `function_map.h`, `nudb_kmer_db.h`, `fasta_parser.*`, `operators.h`,
  `seed_utils.h`). Its `Makefile` adds `-I../signature_kmers/src`, but same-directory
  includes shadow it, so those binaries compile against pattyfam's own stale copies.
  **The stored kmer record is binary-incompatible** between the trees (pattyfam's
  `KData` carries `otu_index`/`weight`/`n_proteins`; this module's `StoredKmerData`
  dropped them), so their databases cannot be interchanged.
- Do **not** try to unify by deleting pattyfam's duplicate headers — that silently
  switches resolution to this module's incompatible record layout. The record format
  and index widths must be reconciled first.

Full details, evidence, and a task list for finishing the migration are in
[`MIGRATION_NOTES.md`](MIGRATION_NOTES.md).

## Conventions

- Templates over `K` and over the DB type are the norm; most logic lives in headers
  and `.tcc` include files rather than `.cc`.
- Programs use Boost.Program_options with both flags and positional args; `data-dir`
  is conventionally the first positional and holds `kmer_data.*` + `function.index`.
- Callbacks (`hit_cb`, `call_cb`) are passed as template functor params, not
  `std::function`, on the hot paths.
