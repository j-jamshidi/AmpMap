# Makefile for ONT Amplicon Phase Pipeline

.PHONY: help install install-dev test lint format clean build docker-build docker-run

# Default target
help:
	@echo "Available targets:"
	@echo "  install      - Install the package"
	@echo "  install-dev  - Install in development mode with dev dependencies"
	@echo "  test         - Run tests"
	@echo "  lint         - Run linting checks"
	@echo "  format       - Format code with black"
	@echo "  clean        - Clean build artifacts"
	@echo "  build        - Build distribution packages"
	@echo "  docker-build - Build Docker image"
	@echo "  docker-run   - Run Docker container"

# Installation targets
install:
	pip install .

install-dev:
	pip install -e ".[dev]"

# Testing and quality
test:
	pytest tests/ -v --cov=ont_amplicon_phase --cov-report=html

lint:
	flake8 src/ont_amplicon_phase/
	mypy src/ont_amplicon_phase/

format:
	black src/ont_amplicon_phase/
	black tests/

# Build and distribution
clean:
	rm -rf build/
	rm -rf dist/
	rm -rf *.egg-info/
	find . -type d -name __pycache__ -exec rm -rf {} +
	find . -type f -name "*.pyc" -delete

build: clean
	python -m build

# Docker targets
docker-build:
	docker build -t ont-amplicon-phase:latest .

docker-run:
	docker-compose up ont-amplicon-phase

docker-dev:
	docker-compose up -d ont-amplicon-dev
	docker-compose exec ont-amplicon-dev bash

# Development helpers
setup-dev: install-dev
	pre-commit install

validate-config:
	ont-amplicon-phase config

validate-sample-sheet:
	ont-amplicon-phase validate sample_sheet.csv