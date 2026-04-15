# Contributing

This repository uses GitHub Actions for CI and `release-please` for automated
versioning, changelog updates, tags, and GitHub releases.

## Development Flow

1. Open a feature branch from `main`.
2. Keep changes scoped to one concern when possible.
3. Make sure the CI workflow passes.
4. Merge into `main` with **Squash and merge**.

`release-please` reads the commit history on `main`. Squash merges keep that
history predictable and make the release notes match the pull request intent.

## Pull Request Titles

Use Conventional Commit style in the pull request title:

- `feat: add subtree visualization helper`
- `fix: correct dynamic leaf regression fixture`
- `docs: rewrite quick start section`
- `refactor: simplify area attribute refresh`
- `test: cover public core Python example`
- `build: adjust wheel smoke test workflow`
- `ci: add release-please automation`

Recommended types:

- `feat`: user-visible feature, API addition, new example, new notebook flow
- `fix`: bug fix or behavioral correction
- `docs`: documentation-only change
- `refactor`: internal restructuring without intended behavior change
- `test`: test-only change
- `build`: packaging or build-system change
- `ci`: GitHub Actions or automation change
- `perf`: performance improvement
- `chore`: maintenance change that should not produce a release note unless
  there is a clear reason

Add `!` after the type for breaking changes:

- `feat!: simplify DynamicComponentTree constructor`

If a change is breaking, describe the migration impact clearly in the pull
request body.

## Merge Policy

- Prefer **Squash and merge** for regular pull requests.
- Keep the squash commit title equal to the pull request title.
- Avoid rewriting the squash title into a non-Conventional Commit message.
- Avoid merge commits for routine feature work, because they make automated
  release notes noisier.

## Release Flow

After changes land on `main`, `release-please` may open or update a release
pull request. When that release PR is merged, the workflow will:

1. update `CHANGELOG.md`;
2. create the version tag;
3. create the GitHub release;
4. build the Python distributions and attach them to the workflow run.

The package version is derived from Git tags via `setuptools-scm`, so release
tags and built artifacts should stay aligned without editing a hardcoded
version string.

## Versioning Convention

This repository follows SemVer through Conventional Commits and
`release-please`.

- `feat:` means a backward-compatible feature and should produce a **minor**
  version bump.
- `fix:` means a bug fix or behavioral correction and should produce a
  **patch** version bump.
- `type!:` or a `BREAKING CHANGE:` footer means an incompatible API or workflow
  change and should produce a **major** version bump.

Use the other common types with narrower intent:

- `docs:` for user-facing documentation changes. In Python repositories,
  `release-please` treats `docs` as releasable, so use it only when the
  documentation change is worth appearing in release notes.
- `perf:` for measurable performance improvements without API breakage.
- `refactor:`, `test:`, `build:`, `ci:`, and `chore:` for maintenance work that
  should normally not define the public versioning intent on their own.

Practical rule for the team:

- choose `feat` when users gain a new capability;
- choose `fix` when users get the same API with corrected behavior;
- choose `!` or `BREAKING CHANGE` when users must adapt code, scripts, or
  workflows;
- choose `docs`, `ci`, `build`, `test`, `refactor`, or `chore` only when the
  change is not primarily a feature or a fix.

Examples:

- `feat: add dynamic subtree update notebook`
- `fix: restore public wheel import on macOS`
- `docs: clarify PyPI installation and release flow`
- `feat!: rename DynamicComponentTreeAdjustment constructor`

If you need to force a specific version outside the normal rules, use a commit
body footer such as `Release-As: 0.3.0` in the release workflow.

## PyPI Publication

PyPI publication is intentionally manual.

After a release PR is merged and the release workflow finishes, publish from a
local, reviewed checkout:

```sh
git fetch --tags
python -m build
python -m twine check dist/*
python -m twine upload dist/*
```

This keeps the final publication step under direct review while preserving
automated versioning and release notes.

## Before Opening a Pull Request

- Run the relevant local tests.
- Update examples or notebooks when the public workflow changes.
- Update documentation when the public API or release behavior changes.
