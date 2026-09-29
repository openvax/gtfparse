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
