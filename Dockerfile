FROM --platform=linux/amd64 python:3.10-slim

WORKDIR /app

# `freesasa` builds a C extension from source if no matching wheel exists
RUN apt-get update && \
  apt-get install -y --no-install-recommends build-essential && \
  rm -rf /var/lib/apt/lists/*

COPY --from=ghcr.io/astral-sh/uv:latest /uv /uvx /usr/local/bin/

COPY . .
RUN chmod +x src/tools/ANARCI/hmmscan src/tools/voronota/voronota

RUN uv sync

ENTRYPOINT ["uv", "run", "deeprank-ab-predict"]
