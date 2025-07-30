"""File handling utilities."""

import os
import shutil
from pathlib import Path
from typing import List, Optional
import logging

logger = logging.getLogger(__name__)


def setup_output_directory(output_dir: Path, barcode: str, clean: bool = False) -> Path:
    """
    Setup output directory structure for a barcode.
    
    Args:
        output_dir: Base output directory
        barcode: Barcode identifier
        clean: Whether to clean existing directory
    
    Returns:
        Path to the barcode-specific directory
    """
    barcode_dir = output_dir / barcode
    
    if clean and barcode_dir.exists():
        logger.info(f"Cleaning existing directory: {barcode_dir}")
        shutil.rmtree(barcode_dir)
    
    barcode_dir.mkdir(parents=True, exist_ok=True)
    logger.info(f"Created output directory: {barcode_dir}")
    
    return barcode_dir


def cleanup_files(file_paths: List[Path], ignore_errors: bool = True) -> None:
    """
    Clean up temporary files and directories.
    
    Args:
        file_paths: List of file/directory paths to remove
        ignore_errors: Whether to ignore deletion errors
    """
    for path in file_paths:
        try:
            if path.exists():
                if path.is_dir():
                    shutil.rmtree(path)
                    logger.debug(f"Removed directory: {path}")
                else:
                    path.unlink()
                    logger.debug(f"Removed file: {path}")
        except Exception as e:
            if not ignore_errors:
                raise
            logger.warning(f"Failed to remove {path}: {e}")


def ensure_file_exists(file_path: Path, description: str = "File") -> None:
    """
    Ensure a file exists, raise FileNotFoundError if not.
    
    Args:
        file_path: Path to check
        description: Description for error message
    """
    if not file_path.exists():
        raise FileNotFoundError(f"{description} not found: {file_path}")


def copy_template_file(template_path: Path, destination: Path, replacements: Optional[dict] = None) -> None:
    """
    Copy a template file to destination with optional text replacements.
    
    Args:
        template_path: Source template file
        destination: Destination file path
        replacements: Dictionary of text replacements {old: new}
    """
    ensure_file_exists(template_path, "Template file")
    
    with open(template_path, 'r') as f:
        content = f.read()
    
    if replacements:
        for old, new in replacements.items():
            content = content.replace(old, str(new))
    
    with open(destination, 'w') as f:
        f.write(content)
    
    logger.debug(f"Copied template {template_path} to {destination}")


def get_file_size(file_path: Path) -> int:
    """Get file size in bytes."""
    return file_path.stat().st_size if file_path.exists() else 0


def create_symlink(source: Path, link: Path, force: bool = False) -> None:
    """
    Create a symbolic link.
    
    Args:
        source: Source file path
        link: Link path to create
        force: Whether to overwrite existing link
    """
    if link.exists() or link.is_symlink():
        if force:
            link.unlink()
        else:
            raise FileExistsError(f"Link already exists: {link}")
    
    link.symlink_to(source)
    logger.debug(f"Created symlink: {link} -> {source}")