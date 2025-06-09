import logging
import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
import gc
import json # For saving list of columns

from core.pipeline_stage import PipelineStage
from core.utils import Utils # Assuming Utils is in the main/core directory

class PreprocessingPipeline(PipelineStage):
    def _run(self):
        df_raw = self.data_manager.get_data('raw')
        if df_raw is None:
            logging.error("Raw data not loaded. Aborting preprocessing.")
            return

        df = df_raw.copy()
        time_index_for_plots = df.index # Default to sample index
        time_label_for_plots = "Sample Index"

        if self.config.TIMESTAMP_COL in df.columns:
            # Ensure the timestamp column is numeric before attempting conversion
            df = pd.to_numeric(df, errors='coerce') 
            df['datetime'] = Utils.ntp_to_datetime_vectorized(df)
            if not df['datetime'].isnull().all():
                time_index_for_plots = df['datetime']
                time_label_for_plots = "Timestamp"
        elif 'datetime' not in df.columns:
             df['datetime'] = pd.Series(*len(df), index=df.index)
        
        self.data_manager.set_object('time_index_for_plots', time_index_for_plots)
        self.data_manager.set_object('time_label_for_plots', time_label_for_plots)
        
        self.data_manager.set_data('original_with_datetime', df.copy())
        self.data_manager.save_data_to_pickle('original_with_datetime', self.config.DATETIME_COL_ORIGINAL_DF_PATH)

        excluded_cols_base = {self.config.TIMESTAMP_COL, 'datetime', 'Event', 'Description'}
        frame_cols = {col for col in df.columns if "Frame" in col}
        known_non_sensor_cols = {
            'Frame_Pupil', 'Frame_GP', 'PupilFrame', 'GoProFrame', 'NPSample',
            'Center of Mass_Frame', 'Markers_Marker Text'
        }
        excluded_cols = excluded_cols_base | frame_cols | known_non_sensor_cols
        
        potential_sensor_cols = [col for col in df.columns if col not in excluded_cols]
        
        temp_sensor_df = df[potential_sensor_cols].copy()
        for col in potential_sensor_cols:
            temp_sensor_df[col] = pd.to_numeric(temp_sensor_df[col], errors='coerce')

        initial_sensor_cols_df = temp_sensor_df.select_dtypes(include=np.number)
        initial_sensor_cols = initial_sensor_cols_df.columns.tolist()

        if not initial_sensor_cols:
            logging.error("No numerical sensor columns found after attempting conversion! Preprocessing cannot continue meaningfully.")
            # Set empty data and objects to prevent downstream errors if pipeline continues
            self.data_manager.set_data('imputed', pd.DataFrame())
            self.data_manager.set_data('scaled', pd.DataFrame())
            self.data_manager.set_object('final_sensor_cols',)
            self.data_manager.set_object('scaler', None)
            # Save empty files or indicators
            self.data_manager.save_data_to_csv('imputed', self.config.PROCESSED_IMPUTED_DF_PATH)
            self.data_manager.save_data_to_csv('scaled', self.config.PROCESSED_SCALED_DF_PATH)
            self.data_manager.save_object_to_joblib('scaler', self.config.SCALER_OBJECT_PATH)
            self.data_manager.save_object_to_json('final_sensor_cols', self.config.FINAL_SENSOR_COLS_PATH)
            # Save an empty time index file
            pd.DataFrame(columns=['sample_index']).to_csv(self.config.TIME_INDEX_PATH, index=False)
            return

        sensor_data_raw = df[initial_sensor_cols].copy()
        sensor_data_cleaned = sensor_data_raw.dropna(axis=1, how='all')
        final_sensor_cols = sensor_data_cleaned.columns.tolist()

        if not final_sensor_cols:
            logging.error("No sensor columns left after cleaning all-NaNs! Preprocessing cannot continue meaningfully.")
            self.data_manager.set_data('imputed', pd.DataFrame())
            self.data_manager.set_data('scaled', pd.DataFrame())
            self.data_manager.set_object('final_sensor_cols',)
            self.data_manager.set_object('scaler', None)
            self.data_manager.save_data_to_csv('imputed', self.config.PROCESSED_IMPUTED_DF_PATH)
            self.data_manager.save_data_to_csv('scaled', self.config.PROCESSED_SCALED_DF_PATH)
            self.data_manager.save_object_to_joblib('scaler', self.config.SCALER_OBJECT_PATH)
            self.data_manager.save_object_to_json('final_sensor_cols', self.config.FINAL_SENSOR_COLS_PATH)
            pd.DataFrame(columns=['sample_index']).to_csv(self.config.TIME_INDEX_PATH, index=False)
            return

        logging.info(f"Processed {len(final_sensor_cols)} sensor columns after cleaning all-NaNs.")
        self.data_manager.set_object('final_sensor_cols', final_sensor_cols)
        del sensor_data_raw; gc.collect()

        logging.info(f"Imputing missing values in {len(final_sensor_cols)} sensor columns using linear interpolate and ffill/bfill...")
        sensor_data_imputed_df = sensor_data_cleaned.copy()
        for col in final_sensor_cols:
            if sensor_data_imputed_df[col].isnull().any():
                sensor_data_imputed_df[col] = sensor_data_imputed_df[col].interpolate(method='linear', limit_direction='both')
                sensor_data_imputed_df[col] = sensor_data_imputed_df[col].fillna(method='ffill').fillna(method='bfill')
        
        remaining_nans = sensor_data_imputed_df.isnull().sum().sum()
        if remaining_nans > 0:
            logging.warning(f"There are still {remaining_nans} NaN values in sensor data after imputation.")
        else:
            logging.info("Imputation complete. No NaNs remaining in sensor data.")
        
        self.data_manager.set_data('imputed', sensor_data_imputed_df)

        scaler = StandardScaler()
        if not sensor_data_imputed_df.empty and sensor_data_imputed_df.shape[1] > 0:
            sensor_values_scaled = scaler.fit_transform(sensor_data_imputed_df)
            sensor_data_scaled_df = pd.DataFrame(sensor_values_scaled, columns=final_sensor_cols, index=sensor_data_imputed_df.index)
            self.data_manager.set_object('scaler', scaler)
        else:
            logging.warning("Imputed data is empty or has no columns. Skipping scaling.")
            sensor_data_scaled_df = pd.DataFrame(columns=final_sensor_cols, index=sensor_data_imputed_df.index)
            self.data_manager.set_object('scaler', None)

        self.data_manager.set_data('scaled', sensor_data_scaled_df)
        del sensor_data_cleaned; gc.collect()
        if 'sensor_values_scaled' in locals(): del sensor_values_scaled; gc.collect()

        self.data_manager.save_data_to_csv('imputed', self.config.PROCESSED_IMPUTED_DF_PATH)
        self.data_manager.save_data_to_csv('scaled', self.config.PROCESSED_SCALED_DF_PATH)
        self.data_manager.save_object_to_joblib('scaler', self.config.SCALER_OBJECT_PATH)
        self.data_manager.save_object_to_json('final_sensor_cols', self.config.FINAL_SENSOR_COLS_PATH)

        if time_index_for_plots is df['datetime'] and not time_index_for_plots.isnull().all():
            # Ensure sensor_data_imputed_df.index is valid for df.loc
            if sensor_data_imputed_df.index.isin(df.index).all():
                aligned_time_index = df.loc[sensor_data_imputed_df.index, 'datetime']
                aligned_time_index.to_csv(self.config.TIME_INDEX_PATH, header=['datetime'], index=False)
                logging.info(f"Saved aligned datetime index to: {self.config.TIME_INDEX_PATH}")
            else:
                logging.warning("Index mismatch between imputed data and original datetime. Saving numerical index.")
                pd.Series(sensor_data_imputed_df.index).to_csv(self.config.TIME_INDEX_PATH, header=['sample_index'], index=False)
                logging.info(f"Saved numerical sample index to: {self.config.TIME_INDEX_PATH}")
        else:
            pd.Series(sensor_data_imputed_df.index).to_csv(self.config.TIME_INDEX_PATH, header=['sample_index'], index=False)
            logging.info(f"Saved numerical sample index to: {self.config.TIME_INDEX_PATH}")