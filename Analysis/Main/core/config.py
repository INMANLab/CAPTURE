import os

class Config:
    """
    Centralized configuration for the data analysis pipeline.
    """
    def __init__(self):
        # Input
        self.RAW_DATA_FILE_PATH = "/home/iris/research/neuroscience/RW/RW1/Walk1/Merged_example_full.csv"
        self.TIMESTAMP_COL = 'NTP_Timestamp'
        self.KNOWN_DTYPES = {
            # self.TIMESTAMP_COL: 'float64', 
            # 'Markers_Marker Text': 'object'
        }

        # Preprocessing Output Paths
        # Use absolute paths or paths relative to a defined project root
        # For simplicity, keeping existing paths
        project_base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")) # Assuming core is one level down from Main

        self.PREPROCESSED_DATA_DIR = os.path.join(project_base_dir, "explore", "preprocessed_data", "full_oop_structured")
        self.PROCESSED_IMPUTED_DF_PATH = os.path.join(self.PREPROCESSED_DATA_DIR, "sensor_data_imputed.csv")
        self.PROCESSED_SCALED_DF_PATH = os.path.join(self.PREPROCESSED_DATA_DIR, "sensor_data_scaled.csv")
        self.SCALER_OBJECT_PATH = os.path.join(self.PREPROCESSED_DATA_DIR, "standard_scaler.joblib")
        self.FINAL_SENSOR_COLS_PATH = os.path.join(self.PREPROCESSED_DATA_DIR, "final_sensor_cols.json")
        self.TIME_INDEX_PATH = os.path.join(self.PREPROCESSED_DATA_DIR, "time_index.csv")
        self.DATETIME_COL_ORIGINAL_DF_PATH = os.path.join(self.PREPROCESSED_DATA_DIR, "original_df_with_datetime.pkl")

        # Analysis Output Paths
        self.ANALYSIS_BASE_OUTPUT_DIR = os.path.join(project_base_dir, "explore", "graphs", "modular_analysis_oop_structured")
        self.LOG_DIR = os.path.join(project_base_dir, "logs")
        self.LOG_FILE_NAME = "modular_analysis_oop_structured.log"

        self.EDA_PLOT_SUBDIR = "eda_plots"
        self.PCA_PLOT_SUBDIR = "pca_plots"
        self.ACF_PLOT_SUBDIR = "acf_plots"

        # Analysis Parameters
        self.PROBLEM_COLUMN_PCA = 'Center of Mass_CoM vel z'
        self.N_PCA_COMPONENTS = 3
        self.MAX_FEATURES_FOR_FULL_CORR_PLOT = 150
        self.MAX_FEATURES_FOR_FULL_CORR_ANNOTATION = 30
        self.MIN_FEATURES_FOR_FOCUSED_CORR_PLOT = 3
        self.ACF_LAGS = 5000
        self.SPECIFIC_COLUMNS_FOR_ACF = None # Example: ['col1', 'col2'] or None for all
        self.CORRELATION_THRESHOLD_FOR_SUBSET = 0.5
        self.MAX_FEATURES_FOR_FOCUSED_CORR = 25
        self.ANNOTATE_FOCUSED_CORR = True
        
        self.NEURAL_CHANNEL_PATTERNS = ['NeuralPace_Ch1', 'NeuralPace_Ch2', 'NeuralPace_Ch3', 'NeuralPace_Ch4']
        self.CORRELATION_THRESHOLD_NEURAL_FOCUS = 0.4

        # Control Flags
        self.RUN_PREPROCESSING = True
        self.RUN_EDA_PLOTS = True
        self.RUN_PCA_ANALYSIS = True
        self.RUN_ACF_ANALYSIS = True
        self.RUN_DATA_AUGMENTATION = False # Example for a new stage