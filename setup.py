#!/usr/bin/env python3

from setuptools import setup, find_packages
from setuptools.command.install import install
import os
import subprocess
import sys

# Read the contents of README file
this_directory = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()


class PostInstallCommand(install):
    """Post-installation for installation mode."""
    def run(self):
        install.run(self)
        # Setup Clair3 environment after installation
        self.setup_clair3()
    
    def setup_clair3(self):
        """Setup Clair3 conda environment."""
        print("Setting up Clair3 environment...")
        try:
            # Check if conda is available
            subprocess.run(['conda', '--version'], check=True, capture_output=True)
            
            # Create Clair3 environment
            cmd = [
                'conda', 'create', '-c', 'conda-forge', '-c', 'bioconda', 
                '-n', 'clair3', 'python=3.9', 'tensorflow=2.15.0', 
                'whatshap', 'samtools', 'parallel', 'xz', 'zlib', 'bzip2', 
                'automake', 'curl', 'pigz', 'cffi', 'make', 'gcc', '-y'
            ]
            subprocess.run(cmd, check=True)
            
            # Setup Clair3
            home_dir = os.path.expanduser('~')
            clair3_dir = os.path.join(home_dir, 'Clair3')
            
            if not os.path.exists(clair3_dir):
                # Clone Clair3
                subprocess.run(['git', 'clone', 'https://github.com/HKU-BAL/Clair3.git', clair3_dir], check=True)
                
                # Get conda prefix for clair3 environment
                result = subprocess.run(['conda', 'info', '--envs'], capture_output=True, text=True)
                clair3_prefix = None
                for line in result.stdout.split('\n'):
                    if 'clair3' in line:
                        clair3_prefix = line.split()[-1]
                        break
                
                if clair3_prefix:
                    # Build Clair3
                    env = os.environ.copy()
                    env['PREFIX'] = clair3_prefix
                    subprocess.run(['make'], cwd=clair3_dir, env=env, check=True)
                    
                    # Download models
                    models_dir = os.path.join(clair3_dir, 'models')
                    os.makedirs(models_dir, exist_ok=True)
                    
                    model_url = 'http://www.bio8.cs.hku.hk/clair3/clair3_models/clair3_models.tar.gz'
                    subprocess.run(['wget', model_url], cwd=models_dir, check=True)
                    subprocess.run(['tar', '-zxf', 'clair3_models.tar.gz'], cwd=models_dir, check=True)
                    subprocess.run(['rm', 'clair3_models.tar.gz'], cwd=models_dir, check=True)
                    
                    # Install pypy
                    bin_dir = os.path.join(clair3_prefix, 'bin')
                    pypy_url = 'https://downloads.python.org/pypy/pypy3.10-v7.3.19-linux64.tar.bz2'
                    subprocess.run(['wget', pypy_url], cwd=bin_dir, check=True)
                    subprocess.run(['tar', '-jxf', 'pypy3.10-v7.3.19-linux64.tar.bz2'], cwd=bin_dir, check=True)
                    
                    pypy_link = os.path.join(bin_dir, 'pypy3')
                    pypy_target = os.path.join(bin_dir, 'pypy3.10-v7.3.19-linux64', 'bin', 'pypy3')
                    if os.path.exists(pypy_target):
                        os.symlink(pypy_target, pypy_link)
                    
                    subprocess.run(['rm', 'pypy3.10-v7.3.19-linux64.tar.bz2'], cwd=bin_dir, check=True)
            
            print(f"✓ Clair3 environment setup complete!")
            print(f"Set ONT_CLAIR3_PATH={clair3_dir} in your environment")
            
        except subprocess.CalledProcessError as e:
            print(f"Warning: Clair3 setup failed: {e}")
            print("You can manually run: ./setup_clair3.sh")
        except FileNotFoundError:
            print("Warning: conda not found. Please install conda and run: ./setup_clair3.sh")

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
    cmdclass={
        'install': PostInstallCommand,
    },
)