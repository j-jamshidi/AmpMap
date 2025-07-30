"""Configuration management utilities."""

import os
import yaml
from pathlib import Path
from typing import Dict, Any, Optional
import logging

logger = logging.getLogger(__name__)


class ConfigManager:
    """Manages configuration loading and validation."""
    
    def __init__(self, config_path: Optional[Path] = None):
        self.config_path = config_path
        self._config: Dict[str, Any] = {}
        self._load_config()
    
    def _load_config(self) -> None:
        """Load configuration from file or use defaults."""
        # Load default config
        default_config_path = Path(__file__).parent.parent / "config" / "default.yaml"
        with open(default_config_path, 'r') as f:
            self._config = yaml.safe_load(f)
        
        # Override with user config if provided
        if self.config_path and self.config_path.exists():
            with open(self.config_path, 'r') as f:
                user_config = yaml.safe_load(f)
                self._deep_update(self._config, user_config)
        
        # Override with environment variables
        self._load_env_overrides()
    
    def _deep_update(self, base_dict: Dict[str, Any], update_dict: Dict[str, Any]) -> None:
        """Deep update dictionary with another dictionary."""
        for key, value in update_dict.items():
            if key in base_dict and isinstance(base_dict[key], dict) and isinstance(value, dict):
                self._deep_update(base_dict[key], value)
            else:
                base_dict[key] = value
    
    def _load_env_overrides(self) -> None:
        """Load configuration overrides from environment variables."""
        env_mappings = {
            'ONT_REFERENCE_GENOME': ['paths', 'reference_genome'],
            'ONT_CLAIR3_PATH': ['paths', 'clair3_path'],
            'ONT_HAPCUT2_PATH': ['paths', 'hapcut2_path'],
            'ONT_BASE_DATA_DIR': ['paths', 'base_data_dir'],
            'ONT_SCRIPT_DIR': ['paths', 'script_dir'],
            'ONT_SOFTWARE_DIR': ['paths', 'software_dir'],
            'ONT_CLAIR3_MODEL': ['paths', 'clair3_model'],
            'ONT_THREADS': ['variant_calling', 'clair3', 'threads'],
            'ONT_MIN_COVERAGE': ['variant_calling', 'clair3', 'min_coverage'],
            'ONT_LOG_LEVEL': ['logging', 'level'],
            'ONT_S3_BUCKET': ['aws', 's3_bucket'],
            'ONT_S3_PREFIX': ['aws', 's3_prefix'],
            'ONT_UPLOAD_S3': ['output', 'upload_to_s3'],
            'ONT_GENERATE_XML': ['output', 'generate_xml'],
            'ONT_BASESPACE_CONFIG': ['basespace', 'config_name'],
        }
        
        for env_var, config_path in env_mappings.items():
            value = os.getenv(env_var)
            if value:
                # Navigate to the nested config location
                current = self._config
                for key in config_path[:-1]:
                    current = current.setdefault(key, {})
                
                # Convert value to appropriate type
                if config_path[-1] == 'threads' or config_path[-1] == 'min_coverage':
                    value = int(value)
                elif config_path[-1] in ['var_pct_full', 'ref_pct_full', 'var_pct_phasing']:
                    value = float(value)
                elif config_path[-1] in ['enable_phasing', 'use_whatshap', 'remove_intermediate', 'upload_to_s3', 'generate_xml']:
                    value = value.lower() in ('true', '1', 'yes', 'on')
                
                current[config_path[-1]] = value
    
    def get(self, key_path: str, default: Any = None) -> Any:
        """Get configuration value using dot notation."""
        keys = key_path.split('.')
        current = self._config
        
        try:
            for key in keys:
                current = current[key]
            return current
        except (KeyError, TypeError):
            return default
    
    def validate_paths(self) -> bool:
        """Validate that required paths exist."""
        required_paths = [
            'paths.reference_genome',
            'paths.clair3_path',
            'paths.hapcut2_path'
        ]
        
        for path_key in required_paths:
            path_value = self.get(path_key)
            if not path_value:
                logger.error(f"Required path not configured: {path_key}")
                return False
            
            path = Path(path_value)
            if not path.exists():
                logger.error(f"Path does not exist: {path_value}")
                return False
        
        return True
    
    @property
    def config(self) -> Dict[str, Any]:
        """Get the full configuration dictionary."""
        return self._config.copy()


def load_config(config_path: Optional[Path] = None) -> ConfigManager:
    """Load configuration from file or environment."""
    return ConfigManager(config_path)