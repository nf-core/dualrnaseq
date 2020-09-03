FROM nfcore/base:1.10.2  #1.9 #FROM nfcore/base:1.10.2, but errors with this line
LABEL authors="Bozena Mika-Gospodorz, Regan Hayward" \
      description="Docker image containing all software requirements for the nf-core/dualrnaseq pipeline"

# Install the conda environment
COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/nf-core-dualrnaseq-1.0/bin:$PATH

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name nf-core-dualrnaseq-1.0 > nf-core-dualrnaseq-1.0.yml
