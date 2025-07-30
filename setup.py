#!/usr/bin/env python3

from setuptools import setup, find_packages
import os

# Read the contents of README file
this_directory = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name="ont-amplicon-phase",
    version="1.0.0",
    author="Clinical Genomics Team",
    author_email="support@example.com",
    description="Professional pipeline for ONT amplicon sequencing analysis and variant phasing",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/your-org/ont-amplicon-phase",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Healthcare Industry",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Medical Science Apps.",
    ],
    python_requires=">=3.8",
    install_requires=[
        "pysam>=0.19.0",
        "pandas>=1.3.0",
        "numpy>=1.21.0",
        "whatshap>=1.4",
        "click>=8.0.0",
        "pyyaml>=6.0",
        "jsonschema>=4.0.0",
        "pathlib2>=2.3.0",
    ],
    extras_require={
        "dev": [
            "pytest>=6.0",
            "pytest-cov>=2.0",
            "black>=21.0",
            "flake8>=3.9",
            "mypy>=0.910",
        ],
    },
    entry_points={
        "console_scripts": [
            "ont-amplicon-phase=ont_amplicon_phase.cli:main",
        ],
    },
    include_package_data=True,
    package_data={
        "ont_amplicon_phase": ["data/*.vcf", "data/*.xml", "config/*.yaml"],
    },
)