# Development

## Setup

```bash
git clone https://github.com/haddocking/deeprank-ab
cd deeprank-ab
uv sync --extra dev   # or: pip install -e ".[dev]"
```

This installs the project in editable mode plus `pytest`. See [CLAUDE.md](CLAUDE.md) for the
pipeline architecture and things to watch for when editing.

## Tests

```bash
uv run pytest tests/ -v   # or: pytest tests/ -v
```

`tests/test_inference.py` is an end-to-end test (`e2e` marker) — it runs the full pipeline
against `example/test.pdb` and checks the predicted DockQ value against a known-good regression
target, not just that the pipeline doesn't crash.

## Docker

Build and test either variant locally before pushing:

```bash
docker build --platform linux/amd64 -f Dockerfile.cpu -t deeprank-ab:cpu .
docker run --rm --user "$(id -u):$(id -g)" -v "$PWD/example":/data -w /data deeprank-ab:cpu test.pdb
```

`--platform linux/amd64` matters even on Apple Silicon — the vendored `hmmscan`/`voronota`
binaries are x86-64 Linux ELF; an arm64 build would hit "exec format error" on them regardless of
host. `Dockerfile.gpu` additionally needs the
[NVIDIA Container Toolkit](https://docs.nvidia.com/datacenter/cloud-native/container-toolkit/latest/install-guide.html)
and `--gpus all` to actually exercise the GPU path — CI can only build-check it (no GPU runners).

## CI

- `ci.yml` — runs the test suite (Python 3.10, Ubuntu) and build-checks both Dockerfiles (no push).
- `docker-publish.yml` — builds and pushes both image variants to GHCR: `pr-<number>-{cpu,gpu}`
  on pull requests (work-in-progress images, not for general use), `<tag>-{cpu,gpu}` +
  `latest-{cpu,gpu}` on release.
- `publish.yml` — builds and publishes the package to PyPI on release, via OIDC trusted
  publishing (no stored token). Requires the trusted publisher to be registered on PyPI first
  (repo + workflow filename + environment name) before the first release.
