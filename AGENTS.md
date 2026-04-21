## Golden Rules

1. **Never commit to `master`.** Always `git checkout -b <feature-branch>` before editing. Land via PR.
2. **Every PR bumps the version.** Even doc-only PRs — at minimum a patch bump. Update `gtfparse/__init__.py::__version__`.
3. **"Done" means merged AND deployed to PyPI** — never stop at merge. After a PR merges, run `./deploy.sh` from a clean master. Skipping deploy = task not done.
4. **File problems as issues, don't silently work around them.** If you hit a bug here or in a sibling openvax/pirl-unc repo, open a GitHub issue on the correct repo and link it from the PR.
5. **After a PR ships, look for the next block of work.** Read open issues across the relevant openvax repos, group by dependency + urgency. Prefer *foundational* changes that unblock multiple downstream improvements; otherwise chain the smallest independent improvements.

---

## Before Completing Any Task

Before considering any code change complete, you MUST:

1. **Run `./lint.sh`** - Verify linting passes (runs `ruff check`)
2. **Run `./test.sh`** - Verify all tests pass with coverage

Do not tell the user you are "done" or that changes are "complete" until both pass.

## Scripts

- `./lint.sh` - Checks linting with ruff (must pass). **Always use this for linting if it exists.**
- `./test.sh` - Runs pytest with coverage (must pass)
- `./lint-and-test.sh` - Runs lint and test back-to-back
- `./deploy.sh` - Deploys to PyPI (gates on lint.sh and test.sh). **Always use this for deploying if it exists.**
- `./develop.sh` - Installs package in development mode (`pip install -e .`)

## Code Style

- Use ruff for linting
- Configuration is in `pyproject.toml` under `[tool.ruff]`
- Target Python version: 3.9+ (CI matrix runs 3.9 / 3.10 / 3.11)

## Project Shape

- Package: `gtfparse/` (thin library — attribute parsing, GTF reading, missing feature creation)
- Tests: `tests/` (pytest; fixtures under `tests/data/`)
- Core deps: polars, pyarrow, pandas (see `requirements.txt`)
- Version lives in `gtfparse/__init__.py` and is exposed via `setuptools.dynamic` in `pyproject.toml`

---

## Workflow Orchestration

### 1. Upfront Planning
- For ANY non-trivial task (3+ steps or architectural decisions): write a detailed spec before touching code
- If something goes sideways, STOP and re-plan immediately — don't keep pushing
- Use planning/verification steps, not just building
- Write detailed specs upfront to reduce ambiguity

### 2. Self-Improvement Loop
- After ANY correction from the user: update `tasks/lessons.md` with the pattern
- Write rules for yourself that prevent the same mistake
- Ruthlessly iterate on these lessons until mistake rate drops
- Review lessons at session start for relevant project

### 3. Verification Before Done
- Never mark a task complete without proving it works
- Diff behavior between the latest code and your changes when relevant
- Ask yourself: "Would a staff engineer approve this?"
- Run tests, check logs, demonstrate correctness

### 4. Demand Elegance (Balanced)
- For non-trivial changes: pause and ask "is there a more elegant way?"
- If a fix feels hacky: "Knowing everything I know now, implement the elegant solution"
- Skip this for simple, obvious fixes — don't over-engineer
- Challenge your own work before presenting it

### 5. Autonomous Bug Fixing
- When given a bug report: just fix it. Don't ask for hand-holding
- Point at logs, errors, failing tests — then resolve them
- Zero context switching required from the user
- Fix failing unit tests without being told how

---

## Task Management

1. **Plan First**: Write plan to `tasks/todo.md` with checkable items
2. **Verify Plan**: Check in before starting implementation
3. **Track Progress**: Mark items complete as you go
4. **Explain Changes**: High-level summary at each step
5. **Document Results**: Add review section to `tasks/todo.md`
6. **Capture Lessons**: Update `tasks/lessons.md` after corrections

---

## Core Principles

- **Simplicity First**: Make every change as simple as possible. Impact minimal code.
- **No Laziness**: Find root causes. No temporary fixes. Senior developer standards.
- **Minimal Impact**: Changes should only touch what's necessary. Avoid introducing bugs.

## Scientific Domain Knowledge
- **Read the literature**: if some code involves genomic or biological concepts (GTF/GFF formats, feature hierarchies, attribute conventions), feel free to search for specs/review papers before changing code that expresses scientific concepts.
- **Flag inconsistencies**: if code expresses a scientific model that's at odds with your understanding, note that inconsistency and ask for clarification.
