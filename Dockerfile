FROM rocker/tidyverse:latest

# Install system dependencies needed for AWS CLI installation
RUN apt-get update && apt-get install -y \
    unzip \
    curl \
    && rm -rf /var/lib/apt/lists/*

# Install AWS CLI v2 using the official method
RUN curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip" \
    && unzip awscliv2.zip \
    && ./aws/install \
    && rm -rf awscliv2.zip aws/

# Install additional R packages not in tidyverse
RUN R -e "install.packages(c('writexl', 'furrr', 'openxlsx'), repos='https://cran.r-project.org/')"

# Install Bioconductor packages
RUN R -e "if (!require('BiocManager', quietly = TRUE)) install.packages('BiocManager'); BiocManager::install('Biostrings')"

# Verify AWS CLI installation
RUN aws --version

WORKDIR /app
CMD ["/bin/bash"]
