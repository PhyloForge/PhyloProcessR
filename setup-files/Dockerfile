# Base image: Micromamba is a lightweight, fast C++ rewrite of Mamba/Conda
FROM mambaorg/micromamba:1.5.8

# Copy the environment.yml into the container
COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml

# Install the environment directly into the base conda environment to avoid activation issues
RUN micromamba install -y -n base -f /tmp/environment.yml && \
    micromamba clean --all --yes

# Set default path so tools can be called directly
ENV PATH="/opt/conda/bin:${PATH}"

WORKDIR /app
CMD ["/bin/bash"]
