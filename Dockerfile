FROM continuumio/anaconda3

WORKDIR /app

RUN apt-get update && \
  apt-get install -y --no-install-recommends build-essential && \
  rm -rf /var/lib/apt/lists/*

COPY environment-gpu.yml .
RUN conda config --system --set channel_priority strict
RUN conda install -n base -y -c conda-forge mamba
RUN mamba env create -f environment-gpu.yml
RUN mamba install -n deeprank-ab -y -c bioconda -c conda-forge anarci --freeze-installed
RUN conda clean -afy

COPY . .
RUN chmod +x src/tools/ANARCI/hmmscan src/tools/voronota/voronota

ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "deeprank-ab", "python3", "scripts/inference.py"]
