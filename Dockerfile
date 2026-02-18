FROM docker.io/rocker/tidyverse:latest
# --- system deps ---
RUN apt-get update && apt-get install -y \
    curl \
    unzip \
    python3 \
    python3-pip \
    jq \
    seqtk \
    parallel \
    && rm -rf /var/lib/apt/lists/*

# --- AWS CLI v2 ---
RUN curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip" \
    && unzip awscliv2.zip \
    && ./aws/install \
    && rm -rf awscliv2.zip aws/

# --- NCBI datasets CLI (datasets + dataformat) ---
# (NCBI publishes static binaries; pin versions if you want reproducibility)
RUN curl -L "https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/LATEST/linux-amd64/datasets" \
      -o /usr/local/bin/datasets \
    && curl -L "https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/LATEST/linux-amd64/dataformat" \
      -o /usr/local/bin/dataformat \
    && chmod +x /usr/local/bin/datasets /usr/local/bin/dataformat

# --- Nextclade CLI ---
# Build with the recent version of Nextclade
RUN curl -L "https://github.com/nextstrain/nextclade/releases/latest/download/nextclade-x86_64-unknown-linux-gnu" \
      -o /usr/local/bin/nextclade \
    && chmod +x /usr/local/bin/nextclade

# --- R packages ---
RUN R -e "install.packages(c('writexl','furrr','openxlsx'), repos='https://cran.r-project.org/')"
RUN R -e "if (!require('BiocManager', quietly=TRUE)) install.packages('BiocManager'); BiocManager::install('Biostrings')"

# --- sanity checks ---
RUN aws --version && datasets --version && seqtk 2>/dev/null || true && nextclade --version

WORKDIR /app
CMD ["/bin/bash"]