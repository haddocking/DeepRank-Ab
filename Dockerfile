FROM continuumio/anaconda3

WORKDIR /app

# Environment first (cacheable layer, changes rarely)
COPY environment-gpu.yml .
RUN conda config --system --set channel_priority strict && \
    conda install -n base -y -c conda-forge mamba && \
    mamba env create -f environment-gpu.yml && \
    mamba install -n deeprank-ab -y -c bioconda -c conda-forge anarci --freeze-installed && \
    conda clean -afy

COPY . .
RUN chmod +x src/tools/ANARCI/hmmscan src/tools/voronota/voronota

ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "deeprank-ab", "python3", "scripts/inference.py"]
