"""Command-line interface for the ONT Amplicon Phase pipeline."""

import click
import logging
import sys
from pathlib import Path
from typing import Optional

from .core.pipeline import AmpliconPipeline
from .utils.config import load_config
from .utils.validators import validate_sample_sheet


def setup_logging(level: str = "INFO") -> None:
    """Setup logging configuration."""
    logging.basicConfig(
        level=getattr(logging, level.upper()),
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        handlers=[
            logging.StreamHandler(sys.stdout),
            logging.FileHandler('ont_amplicon_phase.log')
        ]
    )


@click.group()
@click.option('--config', '-c', type=click.Path(exists=True, path_type=Path),
              help='Path to configuration file')
@click.option('--log-level', default='INFO', 
              type=click.Choice(['DEBUG', 'INFO', 'WARNING', 'ERROR']),
              help='Logging level')
@click.pass_context
def main(ctx: click.Context, config: Optional[Path], log_level: str) -> None:
    """ONT Amplicon Phase Pipeline - Professional bioinformatics pipeline for amplicon analysis."""
    setup_logging(log_level)
    
    # Load configuration
    config_manager = load_config(config)
    
    # Validate configuration
    if not config_manager.validate_paths():
        click.echo("Configuration validation failed. Please check your configuration.", err=True)
        sys.exit(1)
    
    # Store in context for subcommands
    ctx.ensure_object(dict)
    ctx.obj['config'] = config_manager


@main.command()
@click.argument('run_id', type=str)
@click.option('--input-dir', '-i', type=click.Path(exists=True, path_type=Path),
              help='Input directory containing BAM files and sample sheet')
@click.option('--output-dir', '-o', type=click.Path(path_type=Path),
              help='Output directory for results')
@click.option('--sample-sheet', '-s', type=click.Path(exists=True, path_type=Path),
              help='Path to sample sheet CSV file')
@click.option('--clean', is_flag=True, default=False,
              help='Clean existing output directories')
@click.option('--dry-run', is_flag=True, default=False,
              help='Perform dry run without executing pipeline')
@click.option('--no-upload', is_flag=True, default=False,
              help='Skip S3 upload')
@click.option('--no-xml', is_flag=True, default=False,
              help='Skip XML generation')
@click.pass_context
def run(ctx: click.Context, run_id: str, input_dir: Optional[Path], 
        output_dir: Optional[Path], sample_sheet: Optional[Path],
        clean: bool, dry_run: bool, no_upload: bool, no_xml: bool) -> None:
    """Run the complete amplicon analysis pipeline."""
    
    config_manager = ctx.obj['config']
    logger = logging.getLogger(__name__)
    
    # Set default paths if not provided
    if not input_dir:
        input_dir = Path.cwd() / run_id
    if not output_dir:
        output_dir = input_dir / "results"
    if not sample_sheet:
        sample_sheet = input_dir / "sample_sheet.csv"
    
    # Validate inputs
    if not input_dir.exists():
        click.echo(f"Input directory not found: {input_dir}", err=True)
        sys.exit(1)
    
    if not sample_sheet.exists():
        click.echo(f"Sample sheet not found: {sample_sheet}", err=True)
        sys.exit(1)
    
    # Validate sample sheet
    is_valid, errors, df = validate_sample_sheet(sample_sheet)
    if not is_valid:
        click.echo("Sample sheet validation failed:", err=True)
        for error in errors:
            click.echo(f"  - {error}", err=True)
        sys.exit(1)
    
    logger.info(f"Starting pipeline for run: {run_id}")
    logger.info(f"Input directory: {input_dir}")
    logger.info(f"Output directory: {output_dir}")
    logger.info(f"Sample sheet: {sample_sheet}")
    
    if dry_run:
        click.echo("Dry run mode - pipeline would process the following samples:")
        for _, row in df.iterrows():
            click.echo(f"  - {row['Episode']} (Barcode: {row['Barcode']})")
        return
    
    # Override config based on CLI options
    if no_upload:
        config_manager._config['output']['upload_to_s3'] = False
    if no_xml:
        config_manager._config['output']['generate_xml'] = False
    
    # Initialize and run pipeline
    try:
        pipeline = AmpliconPipeline(config_manager.config)
        results = pipeline.run(
            run_id=run_id,
            input_dir=input_dir,
            output_dir=output_dir,
            sample_sheet=sample_sheet,
            clean=clean
        )
        
        # Report results
        click.echo(f"\nPipeline completed successfully!")
        click.echo(f"Processed {len(results)} samples")
        
        for result in results:
            status = "✓" if result['success'] else "✗"
            click.echo(f"  {status} {result['episode']} ({result['barcode']})")
            if not result['success']:
                click.echo(f"    Error: {result.get('error', 'Unknown error')}")
        
    except Exception as e:
        logger.error(f"Pipeline failed: {e}")
        click.echo(f"Pipeline failed: {e}", err=True)
        sys.exit(1)


@main.command()
@click.argument('sample_sheet', type=click.Path(exists=True, path_type=Path))
def validate(sample_sheet: Path) -> None:
    """Validate sample sheet format and content."""
    
    is_valid, errors, df = validate_sample_sheet(sample_sheet)
    
    if is_valid:
        click.echo("✓ Sample sheet validation passed")
        click.echo(f"Found {len(df)} samples to process")
    else:
        click.echo("✗ Sample sheet validation failed:", err=True)
        for error in errors:
            click.echo(f"  - {error}", err=True)
        sys.exit(1)


@main.command()
@click.pass_context
def config(ctx: click.Context) -> None:
    """Display current configuration."""
    
    config_manager = ctx.obj['config']
    config_dict = config_manager.config
    
    click.echo("Current Configuration:")
    click.echo("=" * 50)
    
    def print_config(d, indent=0):
        for key, value in d.items():
            if isinstance(value, dict):
                click.echo("  " * indent + f"{key}:")
                print_config(value, indent + 1)
            else:
                click.echo("  " * indent + f"{key}: {value}")
    
    print_config(config_dict)


@main.command()
def version() -> None:
    """Display version information."""
    from . import __version__, __author__
    click.echo(f"ONT Amplicon Phase Pipeline v{__version__}")
    click.echo(f"Author: {__author__}")


if __name__ == '__main__':
    main()