import os
import logging
import sys
import pandas as pd

class Utils:
    """
    Utility functions for the pipeline.
    """
    @staticmethod
    def setup_logging(log_dir, log_file_name):
        os.makedirs(log_dir, exist_ok=True)
        log_file_path = os.path.join(log_dir, log_file_name)
        
        # Remove existing handlers to avoid duplicate logs if re-running in same session
        for handler in logging.getLogger().handlers:
            logging.root.removeHandler(handler)
            
        logging.basicConfig(
            level=logging.INFO,
            format="%(asctime)s - %(levelname)s - %(module)s.%(funcName)s:%(lineno)d - %(message)s",
            handlers=logging.FileHandler(log_file_path, mode='w'),
        )
        logging.info(f"Logging configured. Log file at: {log_file_path}")

    @staticmethod
    def create_output_directory(base_path, sub_dir_name=None):
        path = os.path.join(base_path, sub_dir_name) if sub_dir_name else base_path
        os.makedirs(path, exist_ok=True)
        logging.info(f"Output directory ensured: {path}")
        return path

    @staticmethod
    def ntp_to_datetime_vectorized(ntp_seconds_series):
        logging.info("Starting vectorized timestamp conversion using NTP epoch (1900-01-01)...")
        if not isinstance(ntp_seconds_series, pd.Series):
            ntp_seconds_series = pd.Series(ntp_seconds_series)
        if ntp_seconds_series.isnull().all():
            logging.warning("Input NTP seconds series is all NaN.")
            return pd.Series(*len(ntp_seconds_series), index=ntp_seconds_series.index)
        try:
            # Ensure the series is numeric before conversion
            numeric_series = pd.to_numeric(ntp_seconds_series, errors='coerce')
            converted_dates = pd.to_datetime(numeric_series, unit='s', origin='1900-01-01', errors='coerce')
            
            # Compare NaNs in original (after to_numeric) vs converted
            if converted_dates.isnull().sum() > numeric_series.isnull().sum():
                logging.warning("Additional NaTs produced during NTP to datetime conversion.")
            logging.info("Timestamp conversion successful.")
            return converted_dates
        except Exception as e:
            logging.error(f"Error during NTP conversion: {e}", exc_info=True)
            return pd.Series(*len(ntp_seconds_series), index=ntp_seconds_series.index)

    @staticmethod
    def sanitize_filename_component(name_str):
        """Sanitizes a string to be used as part of a filename."""
        if not isinstance(name_str, str):
            name_str = str(name_str)
        return name_str.replace(' ', '_').replace('/', '_').replace('.', '_').replace(':', '_')