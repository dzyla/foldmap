# Releasing foldmap

A release is made by pushing a version tag. Do not publish a release from the GitHub web UI: the workflow creates
the GitHub release itself, last, after every gate has passed (`.github/workflows/release.yml`).

1. Set `foldmap.__version__` (`src/foldmap/__init__.py`) to the new version, commit it, and push to `main`. Wait for
   that commit's tests and determinism runs to pass.
2. Tag that commit and push the tag:

       git tag vX.Y.Z
       git push origin vX.Y.Z

3. The release workflow then runs, in this order. Any failure stops everything after it:
   1. **check-tag**: the tag names the package's version (`scripts/check_release.py`);
   2. **tests** (`tests.yml`) and **determinism** (`determinism.yml`: every committed test structure drawn on Linux,
      macOS and Windows, and the figures compared), both on the tagged commit;
   3. **build**: the sdist and the wheel;
   4. **publish** to PyPI by trusted publishing, when the repository variable `PUBLISH_PYPI` is `true`;
   5. **github-release**: the GitHub release, with the built files attached. Zenodo archives the release and mints
      its DOI only now, so a release that failed a gate is never archived.

If a gate fails, fix the problem on `main`, delete the tag (`git push --delete origin vX.Y.Z`, `git tag -d vX.Y.Z`),
and tag the fixed commit with a new version. No release exists yet at that point, so nothing has to be withdrawn.

## The paper's release

The paper describes one frozen release. After its tag's workflow has passed:

- draw and run the fourth random sample once, in a clean checkout of that tag (`docs/validation/fourth-sample-protocol.md`),
  and commit it;
- regenerate the records `docs/validation/freeze-checklist.md` lists;
- regenerate the paper's numbers, tables and figures with `scripts/paper/run_all.sh TAG HARNESS`, HARNESS the commit
  that adds the fourth sample (`scripts/paper/README.md`, "The tag and the harness commit"; the environment is
  `scripts/paper/requirements.lock`).
