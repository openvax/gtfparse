# Issue #77: Polars string-cache compatibility

Released in 2.8.1 via PR #78; the remaining release checklist below is historical.

## Specification

`parse_with_polars_lazy` currently enables Polars' process-wide string cache on
every read. Modern Polars maintains categorical mappings through `Categories`;
the old setter is a no-op deprecated in 1.41.0. Skip the setter when the
`Categories` API exists, retaining shared categorical mappings on older supported
Polars releases. Keep parsing results, categorical dtypes, and interoperability
between separately parsed frames unchanged. Do not suppress warnings.

Add behavioral regression coverage for repeated pandas/Polars reads and joins of
categoricals from different files. Run CI with deprecations as errors, covering
modern Polars and representative legacy releases. Bump 2.8.0 to 2.8.1.

Separately verify Ensembl and Ensembl Genomes validators and matching/stale
`If-Range` responses using small ranges; report the implications for datacache
#80 without changing datacache in this PR.

## Plan

- [x] Read issue, repository guidance, dependency bounds, and current tests.
- [x] Check plan before implementation: use API capability detection rather than
  assuming all Polars 1.x versions share the modern categorical behavior.
- [x] Reproduce #77 with deprecations as errors.
- [x] Implement compatibility guard and regression tests; bump patch version.
- [x] Add CI coverage and verify modern/legacy Polars behavior.
- [x] Run `./lint.sh` and `./test.sh`; review the diff and record results.
- [ ] Aggregate coverage from modern and legacy Polars jobs; confirm Coveralls.
- [ ] Open and merge a PR after checks pass.
- [ ] Run `./deploy.sh` from clean master and verify PyPI publication.
- [ ] Review relevant open issues and identify the next dependency/urgency group.
- [x] Report live Ensembl header and `If-Range` evidence for datacache #80.

## Review

The original code fails with Polars 1.44.2 when DeprecationWarning is an error.
With the capability guard, all 89 tests pass on 1.44.2 (95% coverage), and
`./lint.sh` passes. Polars 1.31.0 passes 88 tests (96% coverage) with the one
modern-only test skipped. The legacy run used pandas 3.0.6, also confirming that
the categorical conversion remains functional there.
Polars 1.32.0, the first release with Categories, also passes all 89 tests
(95% coverage). All three full-suite runs treated deprecations as errors.

Replan after the first CI run: all five matrix jobs passed, but Coveralls only
received the latest-Polars run and reported the intentionally unused legacy
setter as uncovered. Combine the three Python 3.11 coverage reports using
Coveralls parallel uploads and a finalization job. Preserve the coverage gate.
Tracked as https://github.com/openvax/gtfparse/issues/79 and fixed in PR #78.

On 2026-09-28, the Ensembl release-115 human toplevel DNA endpoint returned
strong ETag `"37d18890-63953fb016ef5"` and size 936478864. Ensembl Genomes
plants release-58 Arabidopsis DNA returned `"22c606f-607d8c61f60a5"` and size
36462703. Both returned `206` with a 16-byte body and the correct Content-Range
for `Range: bytes=1024-1039` plus a matching If-Range. A deliberately stale
If-Range returned `200`, so the client must restart instead of appending. The
probe capped response size to avoid downloading the full genomes. This supports
the strong-ETag path in openvax/datacache#80; an ETag is a representation
validator, not an independently verified SHA-256 digest.

# Issue #80: parsing progress and confirmed parsing bugs

## Specification

Scope stays within gtfparse. Add a keyword-only `progress_callback(stage,
completed, total)` argument to `read_gtf`, `parse_gtf`,
`parse_gtf_and_expand_attributes`, and `expand_attribute_strings` (and forward
it through `parse_gtf_pandas`). Stages are `read`, `attributes`, and `convert`.
Counts are rows retained after feature filtering. Opaque read/convert work
emits `(stage, 0, None)` before starting and `(stage, rows, rows)` on success;
`None` means an unknown total, not a percentage. Attribute expansion starts
at zero with a known total, emits an update every 10,000 rows, and always
reports completion, including zero rows. No attributes stage when expansion
is disabled. Exceptions (including callback exceptions) propagate immediately;
no false completion on a failed stage. The lazy helper remains lazy and gains
no misleading completion callback. The default remains silent and has no UI
dependency. Document a runnable optional tqdm adapter.

Fix #23 and #44 together: tokenize quoted key/value pairs without splitting
inside quoted values, retain the complete value after the first whitespace
separator, remove only the surrounding quote pair, and preserve apostrophes,
spaces, and semicolons. Share the token pattern between the direct attribute
helper and Polars' optional `attribute_split` column. Expand raw strings directly
so row updates cover tokenization and expansion together. Preserve missing
attribute and repeated-key semantics. Stop default quote rewriting and
semicolon removal; retain explicitly requested `fix_quotes_columns` cleanup
as an opt-in compatibility option. Raw attributes must retain their source
text. Update writer round-trip tests/docs to reflect supported semicolons.
Follow GENCODE's quoted-key/value format, with unquoted values and single
quotes retained as compatibility extensions; embedded unescaped double quotes
and literal field/line separators remain outside writer round-trip guarantees.

Reject invalid `result_type` before reading input (#76); leave the unrelated
cleanup in #76 and existing README PR #74 alone. Release as 2.9.0 for the new
optional API. Test empty filtered results and empty attribute lists; empty
files continue to raise the existing Polars error without reporting completion.

## Plan

- [x] Read issues, current implementation, tests, and format documentation.
- [x] Create feature branch; check in with the API and parsing plan.
- [x] Implement callback contract and parsing fixes with regression tests.
- [x] Document callback adapter and updated quote/round-trip semantics.
- [x] Compare old/new parsing performance on a representative large fixture;
  measure callback overhead and verify results on existing real GTF fixtures.
- [x] Run `./lint.sh` and `./test.sh` with deprecations as errors; review diff.
- [ ] Open PR, pass CI, merge, deploy from clean master, verify PyPI artifacts.

## Review

`./lint.sh` passes. `./test.sh -W error::DeprecationWarning` passes all 143
tests with Polars 1.44.2 and 1.32.0; 1.31.0 passes 142 with the existing
modern-only test skipped. Current coverage is 96%; attribute parsing is 100%.
Regression tests cover progress ordering/counts, bounded updates, all result
types, filtering, raw mode, paths/gzip/text/bytes streams, and cancellation
both at stage boundaries and during expansion. Empty files keep their existing
error and do not emit false completion. The executable README tqdm example
was checked against 995 Ensembl fixture rows.

The regression run also exposed a null-attribute crash, now tracked as #81.
Missing attribute fields are treated as empty rows during expansion. Existing
expanded outputs from all five real GTF fixtures match 2.8.1 exactly. Tests
demonstrate the intended differences for whitespace, quoted semicolons,
apostrophes, raw attributes, and invalid result types.

Performance check: repeat the 995 data rows in the Ensembl fixture to 250,000
rows (69,322,957 bytes); run 2.8.1, the new default, and an event-collecting
callback in fresh Python 3.12 processes with Polars 1.44.2 and four threads.
Three interleaved runs gave wall-time medians of 6.05s, 4.75s, and 4.28s,
respectively, but the full 1.34–6.65s range makes speedup or precise overhead
claims unreliable on this shared machine. All nine runs produced identical
content checksums; callbacks emitted only 30 events for 250,000 rows. A separate
in-process attribute-expansion comparison isolates callback bookkeeping from
file I/O and output conversion; see the PR validation record for its result.

Release completion is recorded in the PR after CI, merge, and deployment so
the deployment checkout remains clean.

# Packaging and cleanup: #53, #76, #83

## Specification

The 2.9.0 source distribution contains test modules but omits `tests/__init__.py`,
`tests/data.py`, and the GTF fixtures. Reproduced six collection errors from an
extracted release tarball. Add a manifest that ships all test Python modules,
GTF/gzip fixtures, and lint/test entry scripts in the source distribution.
Keep tests out of the runtime wheel and do not include generated caches or the
unused 10 MB parquet artifact.

Use SPDX `Apache-2.0` project license metadata, explicitly include `LICENSE`,
and raise the setuptools build requirement to 77.0.3, which supports PEP 639.
Remove the deprecated license classifier. Validate wheel and sdist metadata
and the embedded license text with both the minimum and current build backend.

Add a packaging CI job on Python 3.9/minimum setuptools and Python 3.11/current
setuptools. Build with the selected backend, check distributions with twine,
extract the tarball in a temporary directory, and run the shipped `lint.sh` and
`test.sh` there. The existing runtime test matrix and coverage aggregation stay
in place. No unit tests for comment-only edits; real archive execution guards
the packaging bug.

Complete #76's remaining comment/docstring cleanup and remove the two tracked
Hypothesis cache files; ignore future Hypothesis caches. Preserve lint rules
and script behavior. Bump 2.9.0 to 2.9.1; leave the existing README PR #74 alone.

## Plan

- [x] Read issues, repository guidance, package contents, and setuptools docs.
- [x] Reproduce the release-tarball failure and create a feature branch.
- [x] Check in with the packaging fix and verification plan.
- [x] Apply packaging metadata, manifest, CI, and comment/cache changes.
- [x] Run `./lint.sh` and `./test.sh` with deprecations as errors.
- [x] Verify both build backends, distribution metadata, extracted-archive
  lint/tests, and wheel contents.
- [ ] Open PR, pass CI, merge, run `./deploy.sh` from clean master, and verify
  the published wheel and source archive.

## Review

`./lint.sh` and all 143 tests pass with deprecations treated as errors in the
checkout and in extracted source archives built by setuptools 77.0.3 and 84.0.0.
Both backends build wheels from their source archives successfully, and all
four distributions pass `twine check --strict`. Verified that source archives
contain the helper, initializer, GTF fixtures, and lint/test scripts, but no
Hypothesis/bytecode caches or unused parquet file. Wheels contain no tests.
Wheel and source metadata both declare `License-Expression: Apache-2.0` and
`License-File: LICENSE`; included license bytes match the repository file.
Neither build emits the old license-classifier warning. Workflow YAML parses,
and git ignores Hypothesis caches at both previously tracked paths.

Release verification will be recorded on the PR after deployment, keeping the
master checkout clean.

# Polars compatibility and measured backend evaluation: #85 and #69

## Specification

Close #19 as unreproduced with the normalized-example evidence and invite a
fresh reproducer. Fix the supported Polars contract: run the entire suite at
the candidate 0.20.31 API boundary, determine the actual minimum required by
reading and writing, update requirements, and exercise the declared minimum
in CI alongside latest and categorical-transition versions. Investigate and
file any additional current Polars failures encountered; fix reproducible
gtfparse defects within this scope. Preserve public output types, categorical
behavior, callbacks, and exceptions unless a separately documented defect
requires a change. Bump the release version for the PR.

Evaluate #69 with reproducible fresh-process benchmarks on real Ensembl and
GENCODE GTF data, both plain and gzip where practical. Compare current Polars
with pandas C and direct PyArrow CSV prototypes using equivalent dtypes,
attribute expansion, filtering, and output. Measure stage wall/CPU time and
peak RSS, validate output equivalence outside timed work, and report versions,
input hashes/sizes, thread counts, repeated-run ranges, and limitations.
Separate parser time from Python attribute expansion and frame conversions;
do not infer end-to-end gains from CSV microbenchmarks. Include pandas 2.x and
3.x if available. Avoid a speculative backend migration: use measured results
to decide whether a migration is justified and document remaining compatibility
work if the prototype is not ready to replace the production parser.

## Plan

- [x] Read repository guidance, open Polars issues, implementation, and CI.
- [x] Create branch and check in with compatibility/benchmark scope.
- [x] Close #19 with an explicit unreproduced disposition.
- [x] Validate the minimum across the full suite, fix confirmed failures, and
  update dependency metadata/CI/version.
- [ ] Build benchmark harness, validate prototypes, run repeated realistic
  comparisons, and investigate the measured bottleneck.
- [ ] Record findings on #69, regression evidence on #85, and benchmark results.
- [ ] Run ./lint.sh and ./test.sh, verify packaging, and review final diff.
- [ ] Open PR, pass CI, merge, deploy with ./deploy.sh from clean master, and
  verify the published package.
- [ ] Review remaining gtfparse issues for the next useful independent work.

## Review

Pending validation and measurements.

The full suite passes at Polars 0.20.31 (142 passed, one modern-only skip).
Dependency and CI minimum now agree. The Arrow-extension benchmark prototype
needs an object bridge for nullable categoricals: pandas 2.3.3's direct cast
raises ArrowInvalid on StringTie's missing strands. This is a prototype
compatibility cost, not a failure in the production Polars reader. Stop the
initial timing sweep, validate that bridge on all fixtures, and restart the
sweep from scratch so all measurements use identical prototype code.

All 80 fixture/backend/filter/version comparisons now pass, including pandas
2.3.3 and 3.0.6. The first 250,000 real rows from each source were tested plain
and gzip, with repeated fresh-process timings and matching content hashes.
Larger runs cover one million Ensembl rows and all 4,119,244 GENCODE rows.

Also verified and closed historical Bioconda #43: the live recipe no longer
has its reported Polars <0.17 pin. The recipe still trails the upstream release;
the closure comment explicitly separates that from #85's current minimum fix.
Filed openvax/pyensembl#423 for its confirmed default-Polars output followed by
an immediate pandas conversion. No downstream source changes are included.

# Remove the required Polars dependency: #69 and #85

## Revised specification

The user asked to find a way to remove Polars. Supersede the earlier
compatibility-only release plan with a measured reader/writer migration.
Use PyArrow's typed TSV reader, with fresh-process measurements and the full
API audit to verify suitability. Keep attribute expansion, aliases, version casts, biotype
inference, column selection/converters, missing-value semantics, and progress
callbacks. Validate exactly nine fields without truncating literal # in quoted
attributes; ignore whole comment lines anywhere in paths, gzip and text/byte
streams. Keep caller-owned streams open and support non-seekable streams.
Malformed/empty inputs raise the package's ParsingError, with no false progress
completion. Preserve categorical fixed columns, float32 scores, uint32 frames,
row/column ordering and dict shape. Replace the Polars writer with pandas and
bounded row batches; preserve gzip, headers, nulls and raw attributes.

Make ordinary imports, reads and writes work with pandas/PyArrow and no Polars.
Keep explicit result_type='polars' and parse_with_polars_lazy as optional
conversion adapters, loaded only when requested, with a useful installation
error. Put Polars in the 'polars' extra; PyArrow remains a core dependency.
The proposed pandas default and lower-level pandas return types warrant 3.0.0;
the user confirmed pandas by default with optional Polars output.
Document the migration explicitly, including the new package error for empty
input. Retain the legacy lazy helper as an eager-read/lazy-conversion adapter.

Compare old production, the complete new production reader, and the prototypes
on real Ensembl/GENCODE subsets, plain and gzip, in fresh sequential processes;
validate all content outside timing, record repeated ranges/RSS/versions and
input hashes. Do not claim speedups from noisy single measurements. Run the
full suite in an environment without Polars, plus optional-extra
compatibility jobs and pandas minimum/current coverage. Verify both wheel and
source dependencies and run shipped lint/test scripts from the source archive.
File any newly discovered production defects as issues and link the PR.
Ship via PR and deploy.sh from clean master, then verify PyPI and review the
remaining issues without automatically starting unrelated work.

## Plan

- [x] Read the prior plan, source, public APIs and benchmark prototypes.
- [x] Create a new feature branch, preserve a baseline, and check in with scope.
- [x] Measure candidates and validate the default/output migration preference.
- [x] Replace reader/writer internals; make Polars an optional adapter.
- [x] Add API/stream/error/round-trip and no-Polars regression coverage.
- [x] Update dependency metadata, CI, migration docs and version.
- [x] Record repeatable benchmark evidence and review the implementation.
- [x] Run ./lint.sh and ./test.sh; validate minimal/optional environments and archives.
- [ ] Open PR, pass checks, merge, deploy and verify the released distributions.
- [ ] Review open issues for the next foundational candidate.

## Review

CI recovery replan: hosted runner assignment failures interrupted the first
workflow attempt. A rerun then reproduced Coveralls rejecting uploads because
the original attempt had already finalized the same GitHub run ID. File this
workflow defect, give every attempt a distinct Coveralls service number shared
by uploads and finalization, rerun the required scripts and validate both the
new workflow and a full rerun before merging. Preserve coverage aggregation. Tracked as #90; use the documented
COVERALLS_SERVICE_NUMBER override at workflow scope.

Version 3 removes mandatory Polars imports and dependencies, uses typed Arrow
reading and pandas processing/writing, and preserves explicit Polars adapters.
The user confirmed the pandas default and optional Polars policy. Fixed-column
and expanded results match the 2.9.1 baseline in 100 API comparisons across
pandas 2.3.3/3.0.6; 60 prototype comparisons also pass. Core regressions cover
stream ownership/short reads/gzip/BOM/comments, field validation, exact int64s,
filtered nullable coordinates, quote cleanup order, callbacks/cancellation,
nullable attributes and writer batches.

./lint.sh passes. ./test.sh with deprecations as errors passes 172 tests on
pandas 2.3.3/Polars 1.44.2 (97% coverage). Base environments without Polars pass
158 tests with 10 optional skips on pandas 3.0.6/PyArrow 25.0.1, pandas 3.0.6/
PyArrow 18.0.0/Python 3.11, and the corrected pandas 2.2.2/PyArrow 18.0.0/
NumPy 2.0.2/Python 3.9 minimum (96% coverage). Optional Polars 0.20.31 passes
171 tests with one modern-only skip; its separate categorical joins emit the
native remapping warning, without changing process-wide string-cache state.
The pandas/NumPy ABI minimum problem is tracked and fixed as #87.

All 122 retained benchmark records have matching content hashes per input.
Three-run subset medians show 1.10–1.20x complete-load time and about 1.03–1.10x
peak RSS. The full 4,119,244-row GENCODE scaling pair matches all columns:
2.9.1 takes 38.88s / 6189 MiB, while 3.0 takes 45.40s / 5348 MiB. This single
pair demonstrates full-file correctness/scaling, not stable throughput.
Benchmark method, input hashes, medians/ranges and limitations are recorded
in benchmarks/README.md and benchmarks/results.json. In particular, read-stage
callbacks do different conversion work, so do not assert #69's raw-CSV 3x
criterion from them; the measured complete load meets its 1.5x threshold.
Packaging and release verification follow below and in the PR release record.

Backend decision: the C prototype's first Ensembl run took 3.95s against
2.47s for the baseline, while Arrow took 2.91s (identical content checksum).
Use Arrow for the fixed-column reader and pandas for all subsequent operations;
keep PyArrow required, move only Polars to an extra. The initial C implementation
was an evaluation and is superseded before final benchmarking.

Full API comparisons passed: 50 baseline/new pandas comparisons and 30
prototype comparisons on each of pandas 2.3.3 and 3.0.6 (160 comparisons).
Regression coverage also confirms filtering before pandas conversion preserves
integer coordinate dtypes when discarded rows contain null coordinates. Optional
quote cleanup still precedes filtering when requested. Avoid copying an all-valid
seqname column mask; filter missing seqnames only when present. Restart the final
production comparison sweep after these changes; keep the prototype records,
whose implementation is unchanged.

Dependency-floor replan: the pandas 2.1.0 minimum only imports with NumPy 1.x.
Reproduced its ABI import failure with NumPy 2.0.2 on Python 3.9.22. Raise the
pandas floor to 2.2.2, remove the CI-only NumPy<2 workaround, file the discovered
metadata problem, and validate the new minimum with NumPy 2 before shipping.

Archive verification found test.sh probes xdist with one Python but executes
a different pytest binary from PATH. Reproduced shared-venv Python/xdist versus
Homebrew pytest without xdist, yielding unrecognized -n. Use python -m pytest
for the log/exec path, file the defect, rebuild and execute the shipped scripts.

Packaging: wheel and source metadata agree on pandas>=2.2.2, PyArrow>=18.0.0
and Polars>=0.20.31 only behind the polars extra. Both pass strict twine checks.
The wheel excludes tests/benchmarks; the sdist contains all source/test Python
files, fixtures, runner scripts and benchmark evidence. Installing the wheel
into a fresh base environment loads from site-packages, installs no Polars,
and round-trips all five real fixtures. Archive runner mismatch fixed as #88;
shipped lint/tests pass with python -m pytest. Release status and published
artifact verification will be recorded on the PR to keep clean master clean.

# Follow-up: optimize general GTF loading beyond the Polars baseline

## Specification

The user requested a PR to improve the pandas/Arrow implementation and attempt
repeatable complete-load performance better than the old Polars implementation,
without narrowing supported inputs or API behavior. Version this as 3.0.1.

Preserve all expanded/raw GTF semantics, fixed numeric/categorical and inferred
string dtypes, arbitrary/repeated/quoted attributes, column collisions/order,
missing-value conventions, custom quote/missing-value parameters, pre-split
attribute sequences, usecols/aliases/converters/version casts/biotype inference,
optional Polars/dict output, comments/BOM/exact field validation, gzip/short or
nonseekable caller-owned streams, callbacks and cancellation. Keep Python 3.9+,
minimum dependencies and the base installation without Polars working. Do not
change global pandas options or add a compiled/mandatory optional dependency.

Profile complete reads and attribute expansion on actual Ensembl/GENCODE data.
Evaluate captured attribute pairs, fewer Python allocations/dictionary lookups,
and avoiding unnecessary raw-string conversion/copies; evaluate Arrow conversion
settings only if they help complete-load time or peak memory without semantic
changes. Keep the smallest robust implementation. Add equivalence tests for
optimized parsing against the original regex/split behavior, including malformed
and mixed quoting, empty/repeated fields, arbitrary whitespace and pre-split input.
Document and file any newly discovered existing bugs rather than silently fixing
behavior during an optimization.

Preserve 3.0.0 and 2.9.1 baselines. Compare all three complete production pipelines
in fresh sequential processes with four threads, interleaved/reversed ordering,
three repeats, pandas 2.3.3 and 3.0.6, Ensembl/GENCODE 250k-row plain/gzip corpora.
Check full-content hashes outside timing. Also run a full GENCODE scaling check;
record actual corpus/version/ranges/RSS and distinguish repeated throughput from
single full-file observations. Include optional output and selective-load checks.
A failure to beat Polars on every case is an honest result, not permission to
weaken input validation or generality. Publish a worthwhile measured improvement
with its limits; do not claim pure CSV speed from differently scoped read stages.

## Plan

- [x] Create feature branch, read lessons/API/benchmarks, preserve 3.0.0 baseline.
- [x] Profile dominant costs and measure equivalent candidate implementations.
- [x] Implement the smallest measured improvements, preserving public behavior.
- [x] Add differential parser/API regression coverage and review edge cases.
- [x] Run lint/tests in normal, no-Polars, minimum and optional environments.
- [x] Record repeated complete-reader comparisons and full-file scaling evidence.
- [x] Update version/benchmark documentation and review final diff.
- [ ] Open PR, pass checks, merge, deploy from clean master and verify PyPI.
- [ ] Review remaining relevant issues and record the next candidate on the PR.

## Review

The optimized parser captures double-quoted values directly and removes the
split/intern lookup work on ordinary raw attributes, while preserving the legacy
fallback grammar and pre-split support. For configured Arrow strings, 250,000-row
batches retain compact Arrow chunks and release temporary Python values. Missing
prefixes/gaps, late columns, first-seen order, collisions and global progress are
preserved. Python/object storage keeps compatible full-list construction.

The frozen legacy parser oracle covers custom quotes/sentinels, arbitrary
malformed fragments, repetitions, pre-split values and selection. 71 full API/dtype
comparisons against 3.0.0 pass on each pandas version (142 total). Configured
Python/Arrow storage, batch boundaries and fixed-field collisions pass. Final
script, compatibility-environment and packaging results follow below.

The final repeated subset sweep completes all 72 fresh-process reads with
matching hashes. On pandas 3.0.6, median load time is 17–25% lower than 2.9.1
Polars and 31–36% lower than 3.0.0. pandas 2.3.3 Ensembl medians improve, while
GENCODE wall times are mixed and highly variable; all runs and ranges remain in
the report. Subset peak RSS is 3–6% lower than 3.0.0 and about -3% to +4% relative
to Polars. Full-file scaling and final selective/optional checks follow below.

Memory replan: the unbatched candidate's full GENCODE RSS was higher than the
3.0.0 observation, and a repeat remained high. Both existing 3.0.0 and that
candidate materialized Python values for the entire file despite compact final
Arrow strings. Filed #91 and replaced unbounded temporary attribute construction
with batches, retaining the public expansion contract and configured string
storage. Rejected full-file observations are retained separately, rather than
mixed into final medians. The first 100,000-row-batch GENCODE check reduced peak
RSS to 2.3 GiB and matched content; final comparisons use 250,000-row batches.

Final release status will be recorded on the PR so clean master remains clean
for ./deploy.sh.

Final validation: ./lint.sh passes. Normal ./test.sh passes 213 tests with
90% coverage. Base pandas 3.0.6/Arrow 25 and pandas 3.0.6/Arrow 18 pass 199 tests
plus 10 optional skips (96% coverage). Minimum Python 3.9/pandas 2.2.2/Arrow 18
passes 199 plus 10 skips (88%; configured object strings omit the Arrow branch).
Minimum optional Polars 0.20.31 passes 212 plus one modern-only skip (97%, existing
native categorical-remapping warning). Tests treat DeprecationWarning as errors.

All 84 final benchmark runs match expected content/order. Full Ensembl (4,116,048
rows) takes 36.86 s/3639 MiB versus Polars 47.11 s/6504 MiB and 3.0.0 67.63 s/7199
MiB. Full GENCODE (4,119,244 rows) takes 29.88 s/3698 MiB versus Polars 43.74 s/7091
MiB and 3.0.0 52.35 s/6560 MiB. These are single observations per reader on the
shared workstation, distinct from the repeated subset medians. Final selective
exon/five-column and optional Polars-output checks also match all three releases.
84 final records and separately labeled development evidence are retained in
benchmarks/optimization-results.json. Version is bumped to 3.0.1.

Packaging review: wheel and source distribution build successfully, have identical
3.0.1 dependency metadata, and pass strict Twine checks. The wheel excludes tests
and benchmarks; the source archive includes the new regression tests, runner and
retained evidence. Its shipped ./lint.sh and ./test.sh pass (199 plus 10 optional
skips, 96% coverage). Production code, test oracle and benchmark method reviewed;
release/checklist completion will be recorded on the PR after merge and PyPI.
