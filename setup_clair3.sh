#!/bin/bash
# Setup script for Clair3 environment

set -e

echo "Setting up Clair3 environment..."

# Create and activate Clair3 environment
conda create -c conda-forge -c bioconda -n clair3 python=3.9 tensorflow=2.15.0 whatshap samtools parallel xz zlib bzip2 automake curl pigz cffi make gcc -y

# Activate environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate clair3

# Clone and build Clair3
if [ ! -d "${HOME}/Clair3" ]; then
    cd ${HOME}
    git clone https://github.com/HKU-BAL/Clair3.git
    cd Clair3
    make PREFIX=${CONDA_PREFIX}
    
    # Download models
    mkdir -p models
    cd models
    wget http://www.bio8.cs.hku.hk/clair3/clair3_models/clair3_models.tar.gz
    tar -zxf clair3_models.tar.gz
    rm clair3_models.tar.gz
fi

# Install pypy
cd ${CONDA_PREFIX}/bin
if [ ! -f "pypy3" ]; then
    wget https://downloads.python.org/pypy/pypy3.10-v7.3.19-linux64.tar.bz2
    tar -jxf pypy3.10-v7.3.19-linux64.tar.bz2
    ln -sf pypy3.10-v7.3.19-linux64/bin/pypy3 pypy3
    rm pypy3.10-v7.3.19-linux64.tar.bz2
fi

echo "Clair3 environment setup complete!"
echo "Set ONT_CLAIR3_PATH=${HOME}/Clair3 in your environment"