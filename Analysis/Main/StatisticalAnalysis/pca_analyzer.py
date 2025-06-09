import logging
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
import os
import gc

from core.pipeline_stage import PipelineStage
from core.utils import Utils

class PCAAnalyzer(PipelineStage):
    def _plot_scree(self, explained_variance_ratio, effective_n_comps, output_dir):
        logging.info(f"Explained variance by {effective_n_comps} PCs: {np.sum(explained_variance_ratio)*100:.2f}%")
        plt.figure(figsize=(8,5))
        plt.plot(range(1, effective_n_comps + 1), explained_variance_ratio, marker='o', linestyle='--')
        plt.title('Scree Plot'); plt.xlabel('Principal Component'); plt.ylabel('Explained Variance Ratio')
        plt.xticks(range(1, effective_n_comps + 1)); plt.grid(); plt.tight_layout()
        plt.savefig(os.path.join(output_dir, "pca_scree_plot.png")); plt.close()
        logging.info("Saved PCA scree plot.")

    def _plot_components_over_time(self, pca_df, time_idx, time_label, effective_n_comps, output_dir):
        pca_plot_time_idx, pca_plot_time_label = time_idx, time_label
        if time_idx is None or len(time_idx)!= len(pca_df):
            logging.warning(f"PCA plot time index mismatch or None. Expected {len(pca_df)}, Got {len(time_idx) if time_idx is not None else 'None'}. Using PCA data index.")
            pca_plot_time_idx = pca_df.index
            pca_plot_time_label = "Sample Index"
        
        plt.figure(figsize=(12,6))
        for i in range(effective_n_comps):
            plt.plot(pca_plot_time_idx, pca_df[f'PC{i+1}'], label=f'PC{i+1}')
        plt.xlabel(pca_plot_time_label); plt.ylabel("Principal Component Value")
        plt.title(f'Top Principal Components Over Time')
        plt.legend(); plt.tight_layout()
        plt.savefig(os.path.join(output_dir, "pca_components_over_time.png")); plt.close()
        logging.info("Saved PCA components plot.")

    def _plot_loadings_heatmap(self, pca, current_pca_cols, effective_n_comps, output_dir):
        if hasattr(pca, 'components_') and pca.components_.shape[1] == len(current_pca_cols):
            loadings = pd.DataFrame(pca.components_.T,
                                    columns=[f'PC{i+1}' for i in range(effective_n_comps)],
                                    index=current_pca_cols)
            plt.figure(figsize=(max(10, len(current_pca_cols)*0.2), max(6, effective_n_comps*0.6)))
            annot = len(current_pca_cols) <= self.config.MAX_FEATURES_FOR_FULL_CORR_ANNOTATION and effective_n_comps <= 10
            sns.heatmap(loadings, annot=annot, cmap='viridis', fmt=".2f")
            plt.title('PCA Loadings'); plt.ylabel('Features'); plt.xlabel('Principal Components')
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, "pca_loadings_heatmap.png")); plt.close()
            logging.info("Saved PCA loadings heatmap.")
            del loadings
        else:
            logging.warning(f"PCA Loadings: Mismatch or no components. Comp features: {pca.components_.shape[1] if hasattr(pca,'components_') else 'N/A'}, Expected: {len(current_pca_cols)}")

    def _run(self):
        output_dir = Utils.create_output_directory(self.config.ANALYSIS_BASE_OUTPUT_DIR, self.config.PCA_PLOT_SUBDIR)
        
        # Define what data this stage needs and how to load it if not present
        # required_data = 
        # Simplified loading for now, direct fetch and file load if None
        data_scaled_df = self.data_manager.get_data('scaled')
        time_idx = self.data_manager.get_object('time_index_for_plots')
        time_label = self.data_manager.get_object('time_label_for_plots')

        if data_scaled_df is None:
            logging.warning("Scaled data not found for PCA. Attempting to load from file.")
            if not self.data_manager.load_data_from_csv('scaled', self.config.PROCESSED_SCALED_DF_PATH):
                logging.error("Failed to load scaled data for PCA. Skipping PCA.")
                return
            data_scaled_df = self.data_manager.get_data('scaled')
        
        if time_idx is None or time_label is None:
            logging.warning("Time index/label not found for PCA. Attempting to load time index.")
            try:
                time_idx_df = pd.read_csv(self.config.TIME_INDEX_PATH)
                if 'datetime' in time_idx_df.columns:
                    time_idx = pd.to_datetime(time_idx_df['datetime'], errors='coerce')
                    time_label = "Timestamp"
                    if time_idx.isnull().all(): 
                        time_idx = data_scaled_df.index if data_scaled_df is not None else pd.Index()
                        time_label = "Sample Index"
                else:
                    time_idx = data_scaled_df.index if data_scaled_df is not None else pd.Index()
                    time_label = "Sample Index"
                self.data_manager.set_object('time_index_for_plots', time_idx)
                self.data_manager.set_object('time_label_for_plots', time_label)
            except Exception as e:
                logging.error(f"Could not load time index for PCA: {e}. Using data index.")
                time_idx = data_scaled_df.index if data_scaled_df is not None else pd.Index()
                time_label = "Sample Index"


        if data_scaled_df is None or data_scaled_df.empty:
            logging.warning("Scaled data for PCA is empty or None. Skipping PCA.")
            return

        data_for_pca = data_scaled_df.copy()
        current_pca_cols = list(data_for_pca.columns)

        if self.config.PROBLEM_COLUMN_PCA and self.config.PROBLEM_COLUMN_PCA in data_for_pca.columns:
            logging.info(f"Excluding '{self.config.PROBLEM_COLUMN_PCA}' from PCA data.")
            data_for_pca = data_for_pca.drop(columns=self.config.PROBLEM_COLUMN_PCA, errors='ignore')
            if self.config.PROBLEM_COLUMN_PCA in current_pca_cols:
                current_pca_cols.remove(self.config.PROBLEM_COLUMN_PCA)
        
        if data_for_pca.empty or data_for_pca.shape[1] == 0:
            logging.warning("No data available for PCA after exclusions. Skipping PCA.")
            return

        max_possible_comps = min(data_for_pca.shape, data_for_pca.shape[1])
        effective_n_comps = min(self.config.N_PCA_COMPONENTS, max_possible_comps)

        if effective_n_comps <= 0:
            logging.warning(f"Cannot perform PCA: effective components = {effective_n_comps}. Skipping.")
            return
        if effective_n_comps < self.config.N_PCA_COMPONENTS:
            logging.warning(f"Reducing N_PCA_COMPONENTS to {effective_n_comps}.")

        pca = PCA(n_components=effective_n_comps)
        try:
            principal_components = pca.fit_transform(data_for_pca)
        except Exception as e:
            logging.error(f"Error during PCA: {e}", exc_info=True)
            return

        pca_df = pd.DataFrame(data=principal_components,
                              columns=[f'PC{i+1}' for i in range(effective_n_comps)],
                              index=data_for_pca.index)
        logging.info(f"PCA result shape: {pca_df.shape}")

        self._plot_scree(pca.explained_variance_ratio_, effective_n_comps, output_dir)
        self._plot_components_over_time(pca_df, time_idx, time_label, effective_n_comps, output_dir)
        self._plot_loadings_heatmap(pca, current_pca_cols, effective_n_comps, output_dir)
        
        del pca_df, principal_components; gc.collect()
