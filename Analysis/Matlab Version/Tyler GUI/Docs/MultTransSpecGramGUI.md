# MultTransSpecGramGUI
Important files for getting started: 
- `Analysis/Matlab Version/Tyler GUI/MultTransSpecGramGUI.mlapp`
- `Analysis/Matlab Version/Tyler GUI/RWAnalysis2.m`
- `RW Task Patient Data/RWAnalysis2_MedianNorm_2sec.mat`

## Using GUI to Plot Spectrograms of Walks

1. **Create an instance** of the `RWAnalysis2` class.  
2. **Set the path** to your `.mat` Analysis File:  
   `RWAnalysis2_MedianNorm_2sec.mat`  
3. **Load the multi-patient data:**  
   `GUIhandle.loadMultData;`  
4. **Launch the spectrogram GUI:**  
   `MultTransSpecGramGUI(GUIhandle);`

Example code:
```matlab
GUIhandle = RWAnalysis2;
GUIhandle.AnalysisFile = '~/MATLAB/Data/RW Task Patient Data/RWAnalysis2_MedianNorm_2sec.mat'; %Replace with local file path to .mat file
GUIhandle.loadMultData;

MultTransSpecGramGUI(GUIhandle);
```
5. **Input the desired parameters in the GUI:**  

![GUI interface](https://github.com/INMANLab/CAPTURE/blob/main/Analysis/Matlab%20Version/Tyler%20GUI/Docs/MultTransSpecGUI_interface.png)
### Core inputs
-  **Save Dir / Save Button**  
  Saves all generated figures and raw data as a `.mat` file with variables `MT` and `MTd` to the specified directory. Automatically creates the folder if it doesn't exist.
- **Transition**  
  Selects the behavioral event to center the analysis (e.g., *Walk Beg*, *Doorway*, *Lost End*).  
  ➤ Use `|` between Transition names to include more than one.

- **Description** *(optional)*  
  Further filters by context (e.g., *Street*, *Stairs*, *ElevatorLobby*).  
  ➤ Use `|` between Description names to include more than one.

- **Region**  
  Brain region(s) to include (e.g., *AntHipp*, *LatTemp*, *Ent+Peri*).  
  ➤ Supports predefined names, `All`, or custom selections. Available channels and patients shown in a table below. 
  ➤ Use `|` between Region names to include more than one.

- **Patient**  
  Choose patient data to include (`1–5`, `All`, or `Leave Out 1–5`).  
  ➤ Manual input of list allowed (e.g., `1,3,4`).

- **Walk**  
  Filter by walking segments (e.g., `All Walks`, `Walk 1–8`, `First/Last Walks`, `Stop Walks`, `Go Walks`).  
  ➤ Manual input of list allowed (e.g., `1,2,7`).

- **Velocity** *(optional)*  
  Filters based on walking speed ±5s from the event, allowing selection of velocity-defined subsets (e.g., `VelHigh`, `VelLow`, or specific terciles/deciles of velocity).

- **Difference**  
  Enables comparison between two sets of conditions (e.g., regions, patients).  
  ➤ Activates second input set; disables baseline correction and single-trial plotting.

### Plot settings:
- **Pval**  
  Significance threshold (default: `0.05`).

- **PvalClust**  
  Cluster-level threshold for permutation testing (default: `0.01`).

- **Permutation**  
  Choose how results are visualized:  
  - `standard`: raw difference (default CLim = `-1,1` or `-10,10` if *Difference* is enabled)  
  - `zscore`: z-scored output (default CLim = `-10,10`)

- **Correction**  
  Choose correction method for multiple comparisons:  
  `cluster`, `pixel`, `fdr`, or `fdr+cluster` (forces z-score mode)  

- **Window**  
  Time range (in seconds) around the transition (default: `-10,10`).

- **Baseline**  
  Time range (in seconds) within `Window` used for normalization.  
  ➤ Disabled when **Full Walk** is selected; normalization is performed across the entire walk segment instead.

- **CLim**  
  Color axis limits for spectrogram plots (default: `-10,10`).

- **YLim**  
  Y-axis range for selected overlaid predictor (e.g., *Vel*, *EDA*, *KDE*).  
  ➤ Auto-set by default; can be customized.

- **Freq (Hz)**  
  Frequency range for plotting individual trial power (used with **Plot Trials**).

- **Box1 / Box2**  
  Define time-frequency windows for comparison before vs. after transition (default Box1: `-4,0,5,8`; Box2: `0,4,5,8`)  
  ➤ Format: **(t1, t2, f1, f2)** = start time, end time, start freq, end freq

- **Plot Power**  
  Generates a boxplot comparing spectral power between Box1 and Box2.

- **Plot Fooof**  
  Generates a boxplot comparing FOOOF-derived spectral metrics.

### Table of electrode regions, sorted by patient and channel numbers.

![Elec_Locs](https://github.com/INMANLab/CAPTURE/blob/main/Analysis/Matlab%20Version/Tyler%20GUI/Docs/ElecLocs.png)

## 
