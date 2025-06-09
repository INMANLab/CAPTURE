import logging
import pandas as pd
import matplotlib.pyplot as plt
from statsmodels.graphics.tsaplots import plot_acf as sm_plot_acf
import os
import gc

from core.pipeline_stage import PipelineStage
from core.utils import Utils

class ACFAnalyzer(PipelineStage):
    def _run(self):
        output_dir = Utils.create_output_directory(self.config.ANALYSIS_BASE_OUTPUT_DIR, self.config.ACF_PLOT_SUBDIR)
        
        data_imputed_df = self.data_manager.get_data('imputed')
        final_sensor_cols = self.data_manager.get_object('final_sensor_cols')

        if data_imputed_df is None:
            logging.warning("Imputed data not found for ACF. Attempting to load from file.")
            if not self.data_manager.load_data_from_csv('imputed', self.config.PROCESSED_IMPUTED_DF_PATH):
                logging.error("Failed to load imputed data for ACF. Skipping ACF.")
                return
            data_imputed_df = self.data_manager.get_data('imputed')
        
        if final_sensor_cols is None:
            logging.warning("Final sensor columns not found for ACF. Attempting to load from file.")
            if not self.data_manager.load_object_from_json('final_sensor_cols', self.config.FINAL_SENSOR_COLS_PATH):
                logging.warning("Failed to load final sensor columns. Using columns from imputed data if available.")
                final_sensor_cols = data_imputed_df.columns.tolist() if data_imputed_df is not None else []
            else:
                final_sensor_cols = self.data_manager.get_object('final_sensor_cols')


        if data_imputed_df is None or data_imputed_df.empty:
            logging.warning("Data for ACF is empty or None. Skipping ACF analysis.")
            return
        if not final_sensor_cols:
             logging.warning("No sensor columns defined for ACF. Skipping.")
             return


        columns_to_analyze = self.config.SPECIFIC_COLUMNS_FOR_ACF if self.config.SPECIFIC_COLUMNS_FOR_ACF else final_sensor_cols
        
        valid_cols_for_acf = [col for col in columns_to_analyze if col in data_imputed_df.columns]
        if not valid_cols_for_acf:
            logging.warning("ACF: No valid columns found after filtering. Skipping.")
            return

        logging.info(f"ACF: Performing for {len(valid_cols_for_acf)} columns, {self.config.ACF_LAGS} lags.")
        for col_name in valid_cols_for_acf:
            series_to_analyze = data_imputed_df[col_name]
            if series_to_analyze.isnull().any():
                logging.warning(f"ACF: NaNs found in imputed '{col_name}'. Dropping for ACF.")
                series_to_analyze = series_to_analyze.dropna()

            if not series_to_analyze.empty and len(series_to_analyze) > 1:
                effective_lags = min(self.config.ACF_LAGS, len(series_to_analyze) - 1)
                if effective_lags <=0:
                    logging.warning(f"Skipping ACF for {col_name}: not enough data points for any lags.")
                    continue

                fig, ax = plt.subplots(figsize=(12, 6))
                try:
                    sm_plot_acf(series_to_analyze, lags=effective_lags, ax=ax, alpha=0.05, use_vlines=True,
                                title=f'Autocorrelation Function (ACF) for {col_name}', fft=True) # fft=True can be faster
                    ax.set_xlabel("Lag"); ax.set_ylabel("Autocorrelation"); plt.tight_layout()
                    plot_filename = f"acf_{Utils.sanitize_filename_component(col_name)}.png"
                    plt.savefig(os.path.join(output_dir, plot_filename)); plt.close(fig)
                    logging.info(f"Saved ACF plot for {col_name}.")
                except Exception as e:
                    logging.error(f"Error plotting ACF for {col_name}: {e}"); plt.close(fig)
            else:
                logging.warning(f"Skipping ACF for {col_name}: insufficient data after NaN handling.")
        gc.collect()