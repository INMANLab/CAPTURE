import os
import logging
import pandas as pd
import joblib
import json

class DataManager:
    """
    Manages data loading, storage, and access throughout the pipeline.
    """
    def __init__(self, config):
        self.config = config
        self.data_store = {} # Stores various dataframes (raw, imputed, scaled, etc.)
        self.objects_store = {} # Stores other objects (scaler, column lists)
        # Ensure preprocessed directory exists if config points to it.
        if hasattr(config, 'PREPROCESSED_DATA_DIR'):
             os.makedirs(self.config.PREPROCESSED_DATA_DIR, exist_ok=True)


    def load_raw_data(self):
        logging.info(f"Loading raw data from {self.config.RAW_DATA_FILE_PATH}...")
        try:
            df = pd.read_csv(self.config.RAW_DATA_FILE_PATH,
                             low_memory=False,
                             dtype=self.config.KNOWN_DTYPES)
            self.data_store['raw'] = df
            logging.info(f"Raw data loaded. Shape: {df.shape}")
            return True
        except FileNotFoundError:
            logging.error(f"Raw data file not found: {self.config.RAW_DATA_FILE_PATH}")
            return False
        except Exception as e:
            logging.error(f"Error loading raw data: {e}", exc_info=True)
            return False

    def get_data(self, key):
        return self.data_store.get(key)

    def set_data(self, key, df):
        self.data_store[key] = df

    def get_object(self, key):
        return self.objects_store.get(key)

    def set_object(self, key, obj):
        self.objects_store[key] = obj

    def save_data_to_csv(self, key, file_path):
        df = self.get_data(key)
        if df is not None:
            try:
                os.makedirs(os.path.dirname(file_path), exist_ok=True)
                df.to_csv(file_path)
                logging.info(f"Saved data '{key}' to: {file_path}")
            except Exception as e:
                logging.error(f"Error saving data '{key}' to CSV: {e}", exc_info=True)
        else:
            logging.warning(f"Data key '{key}' not found in data_store for saving.")

    def load_data_from_csv(self, key, file_path):
        try:
            df = pd.read_csv(file_path, index_col=0) # Assuming index was saved
            self.set_data(key, df)
            logging.info(f"Loaded data '{key}' from: {file_path}")
            return True
        except FileNotFoundError:
            logging.warning(f"File not found for loading data '{key}': {file_path}")
            return False
        except Exception as e:
            logging.error(f"Error loading data '{key}' from CSV: {e}", exc_info=True)
            return False
            
    def save_data_to_pickle(self, key, file_path):
        df = self.get_data(key)
        if df is not None:
            try:
                os.makedirs(os.path.dirname(file_path), exist_ok=True)
                df.to_pickle(file_path)
                logging.info(f"Saved data '{key}' (pickle) to: {file_path}")
            except Exception as e:
                logging.error(f"Error saving data '{key}' to pickle: {e}", exc_info=True)
        else:
            logging.warning(f"Data key '{key}' not found for saving as pickle.")

    def load_data_from_pickle(self, key, file_path):
        try:
            df = pd.read_pickle(file_path)
            self.set_data(key, df)
            logging.info(f"Loaded data '{key}' (pickle) from: {file_path}")
            return True
        except FileNotFoundError:
            logging.warning(f"Pickle file not found for loading data '{key}': {file_path}")
            return False
        except Exception as e:
            logging.error(f"Error loading data '{key}' from pickle: {e}", exc_info=True)
            return False

    def save_object_to_joblib(self, key, file_path):
        obj = self.get_object(key)
        if obj is not None:
            try:
                os.makedirs(os.path.dirname(file_path), exist_ok=True)
                joblib.dump(obj, file_path)
                logging.info(f"Saved object '{key}' to: {file_path}")
            except Exception as e:
                logging.error(f"Error saving object '{key}' with joblib: {e}", exc_info=True)
        else:
            logging.warning(f"Object key '{key}' not found for saving with joblib.")
            
    def load_object_from_joblib(self, key, file_path):
        try:
            obj = joblib.load(file_path)
            self.set_object(key, obj)
            logging.info(f"Loaded object '{key}' from: {file_path}")
            return True
        except FileNotFoundError:
            logging.warning(f"Joblib file not found for loading object '{key}': {file_path}")
            return False
        except Exception as e:
            logging.error(f"Error loading object '{key}' from joblib: {e}", exc_info=True)
            return False

    def save_object_to_json(self, key, file_path):
        obj = self.get_object(key)
        if obj is not None:
            try:
                os.makedirs(os.path.dirname(file_path), exist_ok=True)
                with open(file_path, 'w') as f:
                    json.dump(obj, f, indent=4)
                logging.info(f"Saved object '{key}' to JSON: {file_path}")
            except Exception as e:
                logging.error(f"Error saving object '{key}' to JSON: {e}", exc_info=True)
        else:
            logging.warning(f"Object key '{key}' not found for saving to JSON.")

    def load_object_from_json(self, key, file_path):
        try:
            with open(file_path, 'r') as f:
                obj = json.load(f)
            self.set_object(key, obj)
            logging.info(f"Loaded object '{key}' from JSON: {file_path}")
            return True
        except FileNotFoundError:
            logging.warning(f"JSON file not found for loading object '{key}': {file_path}")
            return False
        except Exception as e:
            logging.error(f"Error loading object '{key}' from JSON: {e}", exc_info=True)
            return False