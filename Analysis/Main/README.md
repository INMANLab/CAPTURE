## Neuroscience Data Analysis Pipeline
This directory contains a modular Python-based pipeline for neuroscience data analysis.

# Folder Structure
core/: Contains core components of the pipeline:
config.py: Centralized configuration.
utils.py: Utility functions (logging, directory creation).
data_manager.py: Handles data loading, saving, and in-memory storage.
pipeline_stage.py: Base class for all analysis stages.
pipeline_orchestrator.py: Manages and runs the sequence of stages.
DataPrep/: Modules for data preparation and preprocessing.
preprocessor.py: Contains the PreprocessingPipeline stage.
DataAugmentation/: Modules for data augmentation techniques (currently a placeholder).
augmentor.py: Placeholder DataAugmentor stage.
StatisticalAnalysis/: Modules for statistical analyses.
pca_analyzer.py: Contains the PCAAnalyzer stage.
acf_analyzer.py: Contains the ACFAnalyzer stage.
Visualization/: Modules for data visualization.
eda_plotter.py: Contains the EDAPlotter stage.
run_main_analysis.py: The main script to execute the pipeline.

# Setup
Ensure all dependencies are installed (e.g., pandas, numpy, scikit-learn, statsmodels, matplotlib, seaborn, joblib).bash pip install pandas numpy scikit-learn statsmodels matplotlib seaborn joblib
Modify core/config.py to set appropriate file paths and parameters for your dataset and analysis. Pay special attention to RAW_DATA_FILE_PATH and output directories.

# Running the Pipeline
Execute the main script from the /research/neuroscience/CAPTURE/Analysis/Main/ directory:

Bash

python run_main_analysis.py
Control which stages are run by modifying the RUN_ flags in core/config.py.

# Adding New Stages
Create a new Python file in the relevant subdirectory (e.g., StatisticalAnalysis/my_new_analysis.py).
Define a new class that inherits from core.pipeline_stage.PipelineStage.
Implement the _run(self) method with the logic for your new analysis stage.
Use self.data_manager to get input data and store results.
Use self.config to access configuration parameters.
Use core.utils.Utils for helper functions.
Import your new stage class in run_main_analysis.py.
Add a control flag for your new stage in core/config.py.
Conditionally add an instance of your new stage to the pipeline in run_main_analysis.py.