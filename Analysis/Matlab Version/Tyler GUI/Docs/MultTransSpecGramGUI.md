Important files for getting started: 
- Analysis > Matlab Version > Tyler GUI >  MultTransSpecramGUI.mlapp
- Analysis > Matlab Version > Tyler GUI > RWAnalysis2.m
- RW Task Patient Data > RWAnalysis2_MedianNorm_2sec.mat 

## Using GUI to Plot Spectrograms of Walks

(1) Create an instance of the `RWAnalysis2` class.
(2) Set the path to your `.mat` Analysis File:  `RWAnalysis2_MedianNorm_2sec.mat`
(3) Load the multi-patient data:  `GUIhandle.loadMultData;`
(4) Launch the spectrogram GUI:  `MultTransSpecGramGUI(GUIhandle);`

```
GUIhandle = RWAnalysis2;
GUIhandle.AnalysisFile = '~/MATLAB/Data/RW Task Patient Data/RWAnalysis2_MedianNorm_2sec.mat'; %Replace with local file path to .mat file
GUIhandle.loadMultData;

MultTransSpecGramGUI(GUIhandle);
```

(5) Input the desired parameters in the GUI: 

![[MultTransSpecGUI_interface.png]]

- **Save Dir / Save Button:** Saves all generated figures and raw data as a .mat file with variables 'MT' and 'MTd' to the specified directory; creates the folder if it doesn't exist.
- **Transition:** Selects the specific behavioral transition during the walk (e.g., Walk Beg/End, Doorways, Lost Beg/End) to center the analysis window and plot trial-aligned data around that event. *| symbol can be added between Transition names to include more than one.* 
- **Description (opt):** (Optional) Further refines the selected transition by context or location (e.g., Street, Stairs, ElevatorLobby). *| symbol can be added between Description names to include more than one.*
- **Region:** Specifies the brain region of interest for analysis (e.g., AntHipp, LatTemp, Ent+Peri, PostHipp+Para). Options include predefined regions, all channels, or a custom selection; available channels and patients shown in a table below. *| symbol can be added between Region names to include more than one.*
- **Patient:** Selects which patient’s data to include in the analysis (Patient 1–5, All Patients, or Leave Out 1–5). *A list of patient numbers can be added manually (e.g. 1,3,4)*
- **Walk:** Filters the data by specific walking segments (e.g., All Walks, Walk 1–8, First/Last Walks, Stop Walks, Go Walks). *A list of walk numbers can be added manually (e.g. 1,2,7)*
- **Velocity (opt):** Filters transitions based on walking speed within ±5 seconds of the event, allowing selection of velocity-defined subsets (e.g., VelHigh, VelLow, or specific terciles/deciles of velocity).
- **Difference:** Enables comparison between two sets of conditions (e.g., transitions, regions, patients) by activating a second set of selection inputs; disables baseline correction and single-trial plotting.
- **Plot settings:**
	- **Pval:** Sets p-value threshold (default = `0.05`). 
	- **PvalClust:** Sets the cluster-level significance threshold (default = `0.01`) for permutation testing.
	- **Permutation:** Selects how permutation test results are visualized—either as raw differences (`standard`, default CLim = `-1,1` or `-10,10` if _Difference_ is checked) or as z-scored values (`zscore`, default CLim = `-10,10`).
	- **Correction:** Chooses the multiple comparisons correction method—`cluster`, `pixel`, `fdr`, or `fdr+cluster` (which forces z-score mode); controls how statistical significance is determined across the spectrogram.
	- **Window:** Sets the time range (in seconds) around the transition event for analysis (default = `-10,10`).
	- **Baseline:** Defines the time window (within the main transition **Window**) used to normalize the spectrogram. If **Full Walk** is selected, this input is disabled and normalization is performed across the entire walk segment instead.
	- **CLim:** Sets the color axis limits (default = `-10,10`) for visualizing z-scored power values in the spectrogram.
	- **YLim:** Sets the y-axis range for the selected overlaid predictor (e.g., Vel, Eda, KDE); default limits auto-adjust by predictor type but can be manually overridden.
	- **Freq (Hz):** Sets the frequency range (in Hz) used for plotting average power of individual trials, which is done only when **Plot Trials** is checked.
	- **Box1**/ **Box2:** Define time-frequency regions (in sec/Hz) for comparing average power before vs. after the transition (defaults: Box1 = `-4,0,5,8`, Box2 = `0,4,5,8`). Each box is specified as **(t1, t2, f1, f2)** = start time, end time, start frequency, end frequency.
	- **Plot Power:** Generates a boxplot comparing average spectral power between Box1 and Box2.
	- **Plot Fooof:** Produces a boxplot comparing FOOOF-derived spectral features (e.g., aperiodic-adjusted power or peak metrics) between Box1 and Box2


![[ElecLocs.png]]

## 