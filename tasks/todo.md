# Issue #77: Polars string-cache compatibility

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

On 2026-09-28, the Ensembl release-115 human toplevel DNA endpoint returned
strong ETag `"37d18890-63953fb016ef5"` and size 936478864. Ensembl Genomes
plants release-58 Arabidopsis DNA returned `"22c606f-607d8c61f60a5"` and size
36462703. Both returned `206` with a 16-byte body and the correct Content-Range
for `Range: bytes=1024-1039` plus a matching If-Range. A deliberately stale
If-Range returned `200`, so the client must restart instead of appending. The
probe capped response size to avoid downloading the full genomes. This supports
the strong-ETag path in openvax/datacache#80; an ETag is a representation
validator, not an independently verified SHA-256 digest.
