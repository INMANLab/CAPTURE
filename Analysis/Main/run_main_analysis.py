import sys
import os
import logging

# Add the project root to the Python path to allow imports from core, etc.
# This assumes run_main_analysis.py is in the 'Main' directory.
PROJECT_ROOT = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, PROJECT_ROOT)

from core.config import Config
from core.utils import Utils
from core.data_manager import DataManager
from core.pipeline_orchestrator import Pipeline

# Import stage classes from their respective modules
from DataPrep.preprocessor import PreprocessingPipeline
# from DataAugmentation.augmentor import DataAugmentor 
from StatisticalAnalysis.pca_analyzer import PCAAnalyzer
from StatisticalAnalysis.acf_analyzer import ACFAnalyzer
from Visualization.eda_plotter import EDAPlotter


def main():
    config = Config()
    Utils.setup_logging(config.LOG_DIR, config.LOG_FILE_NAME) # Pass dir and name separately
    
    data_manager = DataManager(config)
    pipeline = Pipeline(data_manager, config)

    # --- 1. Load Raw Data ---
    if not data_manager.load_raw_data():
        logging.critical("Failed to load raw data. Exiting.")
        sys.exit(1)

    # --- 2. Preprocessing ---
    run_preprocessing_flag = config.RUN_PREPROCESSING
    if not run_preprocessing_flag:
        # Check if all essential preprocessed files exist
        essential_files =  [
            config.PROCESSED_IMPUTED_DF_PATH,
            config.PROCESSED_SCALED_DF_PATH,
            config.SCALER_OBJECT_PATH,
            config.FINAL_SENSOR_COLS_PATH,
            config.TIME_INDEX_PATH,
            config.DATETIME_COL_ORIGINAL_DF_PATH
        ]
        files_exist = all(os.path.exists(f) for f in essential_files)
        if not files_exist:
            logging.info("One or more essential preprocessed files missing and RUN_PREPROCESSING is False. Forcing preprocessing.")
            run_preprocessing_flag = True
        else:
            logging.info("RUN_PREPROCESSING is False and essential preprocessed files exist. Skipping preprocessing stage execution.")
            logging.info("Data will be loaded by analysis stages if needed.")


    if run_preprocessing_flag:
        pipeline.add_stage(PreprocessingPipeline)
    
    # --- Optional: Data Augmentation ---
    if hasattr(config, 'RUN_DATA_AUGMENTATION') and config.RUN_DATA_AUGMENTATION:
        pipeline.add_stage(DataAugmentor) # Assuming DataAugmentor is defined

    # --- 3. EDA Plots ---
    if config.RUN_EDA_PLOTS:
        pipeline.add_stage(EDAPlotter)

    # --- 4. PCA Analysis ---
    if config.RUN_PCA_ANALYSIS:
        pipeline.add_stage(PCAAnalyzer)

    # --- 5. Autocorrelation Analysis ---
    if config.RUN_ACF_ANALYSIS:
        pipeline.add_stage(ACFAnalyzer)
    
    pipeline.run()

if __name__ == "__main__":
    main()