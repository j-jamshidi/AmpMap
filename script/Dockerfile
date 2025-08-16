FROM python:3.9-slim

# Install system dependencies
RUN apt-get update && apt-get install -y \
    samtools \
    && rm -rf /var/lib/apt/lists/*

# Install Python dependencies
RUN pip install pysam numpy

# Copy all Python scripts
COPY *.py /app/

# Set working directory
WORKDIR /app

# Make scripts executable
RUN chmod +x *.py

# Default entrypoint
ENTRYPOINT ["python"]