import logging
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import gc

from core.pipeline_stage import PipelineStage
from core.utils import Utils

class EDAPlotter(PipelineStage):
    def _run(self):
        output_dir = Utils.create_output_directory(self.config.ANALYSIS_BASE_OUTPUT_DIR, self.config.EDA_PLOT_SUBDIR)
        
        data_imputed_df = self.data_manager.get_data('imputed')
        data_scaled_df = self.data_manager.get_data('scaled')
        time_idx = self.data_manager.get_object('time_index_for_plots')
        time_label = self.data_manager.get_object('time_label_for_plots')

        # Attempt to load data if not present in memory
        if data_imputed_df is None:
            if not self.data_manager.load_data_from_csv('imputed', self.config.PROCESSED_IMPUTED_DF_PATH):
                logging.error("Failed to load imputed data for EDA. EDA may be incomplete.")
            data_imputed_df = self.data_manager.get_data('imputed') # re-fetch

        if data_scaled_df is None:
            if not self.data_manager.load_data_from_csv('scaled', self.config.PROCESSED_SCALED_DF_PATH):
                logging.error("Failed to load scaled data for EDA. EDA may be incomplete.")
            data_scaled_df = self.data_manager.get_data('scaled') # re-fetch
        
        if time_idx is None or time_label is None:
            logging.warning("Time index/label not found for EDA. Attempting to load time index.")
            try:
                time_idx_df = pd.read_csv(self.config.TIME_INDEX_PATH)
                if 'datetime' in time_idx_df.columns:
                    time_idx = pd.to_datetime(time_idx_df['datetime'], errors='coerce')
                    time_label = "Timestamp"
                    if time_idx.isnull().all(): 
                        time_idx = data_imputed_df.index if data_imputed_df is not None else pd.Index()
                        time_label = "Sample Index"
                else:
                    time_idx = data_imputed_df.index if data_imputed_df is not None else pd.Index()
                    time_label = "Sample Index"
                self.data_manager.set_object('time_index_for_plots', time_idx)
                self.data_manager.set_object('time_label_for_plots', time_label)
            except Exception as e:
                logging.error(f"Could not load time index for EDA: {e}. Using data index.")
                time_idx = data_imputed_df.index if data_imputed_df is not None else pd.Index()
                time_label = "Sample Index"

        # Use columns from the actual scaled dataframe if available, otherwise empty list
        final_sensor_cols = data_scaled_df.columns.tolist() if data_scaled_df is not None and not data_scaled_df.empty else []

        self._generate_summary_stats(data_imputed_df)
        self._plot_distributions(data_scaled_df, final_sensor_cols, output_dir)
        
        full_corr_matrix = None
        if data_scaled_df is not None and not data_scaled_df.empty:
             full_corr_matrix = self._plot_correlation_matrix_full(data_scaled_df, output_dir)
             self._plot_correlation_matrix_focused(data_scaled_df, full_corr_matrix, output_dir)
             self._plot_correlation_matrix_neural(data_scaled_df, full_corr_matrix, output_dir)
        
        self._plot_time_series(data_scaled_df, final_sensor_cols, time_idx, time_label, output_dir)
        if full_corr_matrix is not None: del full_corr_matrix; gc.collect()


    def _generate_summary_stats(self, data_imputed_df):
        logging.info("Summary Statistics (imputed, original scales):")
        if data_imputed_df is not None and not data_imputed_df.empty:
            logging.info(f"\n{data_imputed_df.describe().T.to_string()}")
        else:
            logging.warning("Skipping summary statistics: imputed data is empty or None.")

    def _plot_distributions(self, data_scaled_df, sensor_cols_list, output_dir):
        logging.info("Plotting distributions (scaled)...")
        if data_scaled_df is None or data_scaled_df.empty:
            logging.warning("Skipping distribution plots: scaled data is empty or None.")
            return
        if not sensor_cols_list: # sensor_cols_list could be an empty list now
            logging.warning("No sensor columns provided for distribution plots.")
            return

        cols_to_plot_dist = list(sensor_cols_list)
        if len(cols_to_plot_dist) > 50:
            logging.warning(f"Plotting distributions for a subset of {len(cols_to_plot_dist)} sensor columns.")
            cols_to_plot_dist = cols_to_plot_dist[:20] + cols_to_plot_dist[-10:]
            cols_to_plot_dist = sorted(list(set(cols_to_plot_dist)))

        for col_name in cols_to_plot_dist:
            if col_name not in data_scaled_df.columns:
                logging.warning(f"Column {col_name} not in scaled data for distribution plot. Skipping.")
                continue
            plt.figure(figsize=(10, 4))
            sns.histplot(data_scaled_df[col_name], kde=True, bins=50)
            plt.title(f'Distribution of {col_name} (Scaled)')
            plt.xlabel("Standardized Value"); plt.ylabel("Frequency")
            plot_filename = f"dist_{Utils.sanitize_filename_component(col_name)}.png"
            plt.savefig(os.path.join(output_dir, plot_filename)); plt.close()
        logging.info(f"Finished plotting distributions for {len(cols_to_plot_dist)} columns.")
        gc.collect()

    def _plot_correlation_matrix_full(self, data_scaled_df, output_dir):
        logging.info("Calculating and plotting FULL correlation matrix...")
        if data_scaled_df is None or data_scaled_df.shape[1] < 2:
            logging.warning("Skipping FULL correlation matrix: scaled data empty, None, or <2 columns.")
            return None
        
        num_features = data_scaled_df.shape[1]
        full_corr_matrix = None

        if num_features > self.config.MAX_FEATURES_FOR_FULL_CORR_PLOT:
            logging.warning(f"Skipping FULL correlation matrix plot: {num_features} features > max {self.config.MAX_FEATURES_FOR_FULL_CORR_PLOT}.")
            if num_features <= 2 * self.config.MAX_FEATURES_FOR_FULL_CORR_PLOT: # Heuristic
                 full_corr_matrix = data_scaled_df.corr()
            else:
                 logging.warning("Skipping calculation of extremely large full_corr_matrix too.")
        else:
            full_corr_matrix = data_scaled_df.corr()
            plt.figure(figsize=(max(12, num_features * 0.3), max(10, num_features * 0.25)))
            annot_full_corr = num_features <= self.config.MAX_FEATURES_FOR_FULL_CORR_ANNOTATION
            sns.heatmap(full_corr_matrix, annot=annot_full_corr, cmap='coolwarm', linewidths=.1, fmt=".1f")
            plt.title('Full Correlation Matrix (Scaled Data)'); plt.tight_layout()
            plt.savefig(os.path.join(output_dir, "correlation_matrix_full.png")); plt.close()
            logging.info("Saved FULL correlation matrix plot.")
        return full_corr_matrix
        
    def _plot_correlation_matrix_focused(self, data_scaled_df, full_corr_matrix, output_dir):
        logging.info("Identifying and plotting FOCUSED correlation matrix...")
        if data_scaled_df is None or data_scaled_df.shape[1] < self.config.MIN_FEATURES_FOR_FOCUSED_CORR_PLOT:
            logging.warning("Skipping focused correlation: data empty/None or too few features.")
            return

        if full_corr_matrix is None:
            logging.info("Full correlation matrix not pre-calculated, calculating for focused subset.")
            if data_scaled_df.shape[1] > 2 * self.config.MAX_FEATURES_FOR_FULL_CORR_PLOT:
                logging.warning("Dataset too large to calculate full_corr_matrix for focused plot. Skipping.")
                return
            full_corr_matrix = data_scaled_df.corr()
            if full_corr_matrix is None: 
                logging.warning("Could not compute full_corr_matrix. Skipping focused plot.")
                return

        abs_corr_upper = full_corr_matrix.abs().where(np.triu(np.ones(full_corr_matrix.shape), k=1).astype(bool))
        selected_features_for_focus = [] # Fixed: Initialize as empty list

        highly_corr_pairs_mask = abs_corr_upper > self.config.CORRELATION_THRESHOLD_FOR_SUBSET
        candidate_features_from_thresh = set()
        for feature_name in abs_corr_upper.columns:
            if highly_corr_pairs_mask[feature_name].any():
                candidate_features_from_thresh.add(feature_name)
                candidate_features_from_thresh.update(abs_corr_upper.index[highly_corr_pairs_mask[feature_name]].tolist())
        
        candidate_features_list = sorted(list(candidate_features_from_thresh))

        if len(candidate_features_list) >= self.config.MIN_FEATURES_FOR_FOCUSED_CORR_PLOT:
            if len(candidate_features_list) > self.config.MAX_FEATURES_FOR_FOCUSED_CORR:
                temp_corr_for_pruning = abs_corr_upper.loc[candidate_features_list, candidate_features_list]
                high_corr_counts = temp_corr_for_pruning.gt(self.config.CORRELATION_THRESHOLD_FOR_SUBSET).sum(axis=1)
                selected_features_for_focus = high_corr_counts.nlargest(self.config.MAX_FEATURES_FOR_FOCUSED_CORR).index.tolist()
            else:
                selected_features_for_focus = candidate_features_list
        
        # This block should only execute if the previous block didn't populate selected_features_for_focus sufficiently
        if len(selected_features_for_focus) < self.config.MIN_FEATURES_FOR_FOCUSED_CORR_PLOT:
            logging.info(f"Threshold method yielded {len(selected_features_for_focus)} features. Selecting from top absolute correlations.")
            all_corr_pairs_sorted = abs_corr_upper.stack().sort_values(ascending=False)
            top_n_candidate_features = set()
            if not all_corr_pairs_sorted.empty:
                for (feat1, feat2), corr_val in all_corr_pairs_sorted.items():
                    if corr_val == 0 and len(top_n_candidate_features) >= self.config.MIN_FEATURES_FOR_FOCUSED_CORR_PLOT : break
                    top_n_candidate_features.add(feat1); top_n_candidate_features.add(feat2)
                    if len(top_n_candidate_features) >= self.config.MAX_FEATURES_FOR_FOCUSED_CORR: break
            
            if len(top_n_candidate_features) >= self.config.MIN_FEATURES_FOR_FOCUSED_CORR_PLOT:
                selected_features_for_focus = sorted(list(top_n_candidate_features)) # Ensure max limit
            elif len(top_n_candidate_features) >= 2:
                selected_features_for_focus = sorted(list(top_n_candidate_features))
            # else: selected_features_for_focus remains (or whatever it was if the first block partially filled it but < MIN_FEATURES)

        if len(selected_features_for_focus) >= 2:
            focused_corr_matrix = data_scaled_df[selected_features_for_focus].corr()
            plt.figure(figsize=(max(10, len(selected_features_for_focus)*0.5), max(8, len(selected_features_for_focus)*0.4)))
            sns.heatmap(focused_corr_matrix, annot=self.config.ANNOTATE_FOCUSED_CORR, cmap='coolwarm', linewidths=.5, fmt=".2f", vmin=-1, vmax=1)
            plt.title(f'Focused Correlation Matrix (Max {len(selected_features_for_focus)} Features)')
            plt.tight_layout(); plt.savefig(os.path.join(output_dir, "correlation_matrix_focused.png")); plt.close()
            logging.info("Saved FOCUSED correlation matrix plot.")
        else:
            logging.info("Not enough features selected (<2) for focused correlation matrix.")

    def _plot_correlation_matrix_neural(self, data_scaled_df, full_corr_matrix, output_dir):
        logging.info("Identifying and plotting NEURAL-FOCUSED correlation matrix...")
        if data_scaled_df is None or data_scaled_df.empty:
            logging.warning("Skipping neural-focused correlation: data empty or None.")
            return

        if full_corr_matrix is None:
            if data_scaled_df.shape[1] > 2 * self.config.MAX_FEATURES_FOR_FULL_CORR_PLOT:
                logging.warning("Dataset too large to calculate full_corr_matrix for neural plot. Skipping.")
                return
            full_corr_matrix = data_scaled_df.corr()
            if full_corr_matrix is None:
                logging.warning("Could not compute full_corr_matrix. Skipping neural plot.")
                return

        all_available_columns = data_scaled_df.columns.tolist()
        actual_neural_channels = [] # Fixed: Initialize as empty list
        for pattern in self.config.NEURAL_CHANNEL_PATTERNS:
            for col in all_available_columns:
                if pattern in col and col not in actual_neural_channels:
                    actual_neural_channels.append(col)
        actual_neural_channels = sorted(list(set(actual_neural_channels)))

        if not actual_neural_channels:
            logging.warning("No neural channels found. Skipping neural-focused correlation matrix.")
            return
        
        features_for_neural_plot = set(actual_neural_channels)
        other_features = [col for col in all_available_columns if col not in actual_neural_channels and col in full_corr_matrix.columns]
        selected_features_for_neural_plot = [] # Fixed: Initialize as empty list

        if other_features and actual_neural_channels: 
            valid_neural_channels_for_corr = [nc for nc in actual_neural_channels if nc in full_corr_matrix.columns]
            if not valid_neural_channels_for_corr:
                logging.warning("Neural channels not found in correlation matrix columns. Skipping some neural logic.")
            else:
                # Ensure other_features are also valid for.loc
                valid_other_features = [of for of in other_features if of in full_corr_matrix.index]
                if valid_other_features:
                    corr_other_to_neural = full_corr_matrix.loc[valid_other_features, valid_neural_channels_for_corr].abs()
                    initial_correlates = set()
                    for other_feat in valid_other_features: 
                        if (corr_other_to_neural.loc[other_feat] > self.config.CORRELATION_THRESHOLD_NEURAL_FOCUS).any():
                            initial_correlates.add(other_feat)
                    features_for_neural_plot.update(initial_correlates)
        
        candidate_list_neural = sorted(list(features_for_neural_plot))

        if len(candidate_list_neural) > self.config.MAX_FEATURES_FOR_FOCUSED_CORR:
            selected_features_for_neural_plot = list(actual_neural_channels) # Start with neural channels
            remaining_slots = self.config.MAX_FEATURES_FOR_FOCUSED_CORR - len(selected_features_for_neural_plot)
            if remaining_slots > 0 and 'corr_other_to_neural' in locals() and not corr_other_to_neural.empty:
                # Filter other_candidates to those present in corr_other_to_neural.index
                other_candidates = [feat for feat in candidate_list_neural if feat not in actual_neural_channels and feat in corr_other_to_neural.index]
                if other_candidates:
                    max_corr_values = corr_other_to_neural.loc[other_candidates].max(axis=1).sort_values(ascending=False)
                    selected_features_for_neural_plot.extend(max_corr_values.head(remaining_slots).index.tolist())
        else:
            selected_features_for_neural_plot = candidate_list_neural

        if len(selected_features_for_neural_plot) < self.config.MIN_FEATURES_FOR_FOCUSED_CORR_PLOT and actual_neural_channels:
            logging.info(f"Neural-focused selection has {len(selected_features_for_neural_plot)} features. Expanding...")
            selected_features_for_neural_plot = set(actual_neural_channels) 
            valid_other_features_for_expansion = [of for of in other_features if of in full_corr_matrix.index] # from earlier
            valid_neural_channels_for_expansion = [nc for nc in actual_neural_channels if nc in full_corr_matrix.columns]

            if valid_other_features_for_expansion and valid_neural_channels_for_expansion:
                max_abs_corr_to_any_neural = full_corr_matrix.loc[valid_other_features_for_expansion, valid_neural_channels_for_expansion].abs().max(axis=1)
                top_other_correlates = max_abs_corr_to_any_neural.sort_values(ascending=False)
                for feat in top_other_correlates.index:
                    selected_features_for_neural_plot.add(feat)
                    if len(selected_features_for_neural_plot) >= self.config.MAX_FEATURES_FOR_FOCUSED_CORR: break
                    if len(selected_features_for_neural_plot) >= self.config.MIN_FEATURES_FOR_FOCUSED_CORR_PLOT and \
                       top_other_correlates[feat] < 0.1 : break # Stop if correlations become too weak
            selected_features_for_neural_plot = sorted(list(selected_features_for_neural_plot))


        if len(selected_features_for_neural_plot) >= 2:
            # Ensure all selected features are actually in data_scaled_df before trying to slice
            valid_selected_features = [f for f in selected_features_for_neural_plot if f in data_scaled_df.columns]
            if len(valid_selected_features) < 2:
                logging.info("Not enough valid features remaining for neural-focused correlation matrix after final check.")
                return

            neural_focused_corr_matrix = data_scaled_df[valid_selected_features].corr()
            plt.figure(figsize=(max(10, len(valid_selected_features)*0.55), max(8, len(valid_selected_features)*0.45)))
            sns.heatmap(neural_focused_corr_matrix, annot=self.config.ANNOTATE_FOCUSED_CORR, cmap='coolwarm_r', linewidths=.5, fmt=".2f", vmin=-1, vmax=1)
            plt.title(f'Neural-Focused Correlation Matrix ({len(valid_selected_features)} Features)')
            plt.tight_layout(); plt.savefig(os.path.join(output_dir, "correlation_matrix_neural_focused.png")); plt.close()
            logging.info("Saved NEURAL-FOCUSED correlation matrix plot.")
        else:
            logging.info("Not enough features selected (<2) for neural-focused correlation matrix.")


    def _plot_time_series(self, data_scaled_df, sensor_cols_list, time_idx, time_label, output_dir):
        logging.info("Plotting time series (scaled)...")
        if data_scaled_df is None or data_scaled_df.empty:
            logging.warning("Skipping time series plots: scaled data is empty or None.")
            return
        if not sensor_cols_list: # sensor_cols_list could be an empty list
            logging.warning("No sensor columns for time series plots.")
            return
        if time_idx is None or len(time_idx) == 0: # Check if time_idx is empty
            logging.warning("Time index not available or empty for time series plots. Using data index.")
            time_idx = data_scaled_df.index
            time_label = "Sample Index"

        cols_to_plot_ts = list(sensor_cols_list)
        if len(cols_to_plot_ts) > 50:
            subset_indices = list(range(min(20, len(cols_to_plot_ts)))) + \
                             list(range(max(20, len(cols_to_plot_ts)-10), len(cols_to_plot_ts)))
            cols_to_plot_ts = [sensor_cols_list[i] for i in sorted(list(set(subset_indices)))]
        
        for col_name_ts in cols_to_plot_ts:
            if col_name_ts not in data_scaled_df.columns:
                logging.warning(f"Column {col_name_ts} not in scaled data for time series plot. Skipping.")
                continue
            plt.figure(figsize=(12, 4))
            data_to_plot = data_scaled_df[col_name_ts]
            
            current_time_axis = time_idx
            current_time_label = time_label
            
            if len(time_idx)!= len(data_to_plot):
                logging.warning(f"Time series for '{col_name_ts}': Length mismatch time index ({len(time_idx)}) vs data ({len(data_to_plot)}).")
                if isinstance(time_idx, pd.DatetimeIndex) and isinstance(data_to_plot.index, pd.RangeIndex):
                    # Try to align if data_to_plot.index is a simple range that might correspond to original positions
                    try:
                        current_time_axis = time_idx[data_to_plot.index]
                        logging.info(f"Aligned datetime index for '{col_name_ts}'.")
                    except IndexError:
                        logging.warning(f"Could not align datetime index for '{col_name_ts}'. Using data's numerical index.")
                        current_time_axis = data_to_plot.index
                        current_time_label = "Sample Index (Fallback)"
                else: # Default to data's own index if alignment is unclear or time_idx is not datetime
                    logging.warning(f"Using data's numerical index for '{col_name_ts}'.")
                    current_time_axis = data_to_plot.index
                    current_time_label = "Sample Index (Fallback)"

            plt.plot(current_time_axis, data_to_plot, label=col_name_ts)
            plt.title(f'Time Series of {col_name_ts} (Scaled)')
            plt.ylabel("Standardized Value"); plt.xlabel(current_time_label)
            plt.legend(); plt.tight_layout()
            plot_filename = f"ts_{Utils.sanitize_filename_component(col_name_ts)}.png"
            plt.savefig(os.path.join(output_dir, plot_filename)); plt.close()
        logging.info(f"Finished plotting {len(cols_to_plot_ts)} time series examples.")
        gc.collect()