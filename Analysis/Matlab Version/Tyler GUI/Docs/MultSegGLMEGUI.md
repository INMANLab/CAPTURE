# MultSegGLMEGUI
Description: Fits a Generalized Linear Mixed-Effects (GLME) model to segments of walking trial data. 

**NOTE:** To use this, make sure you have the Statistics and Machine Learning Toolbox installed in MATLAB. 

Important files for getting started: 
- `Analysis/Matlab Version/Tyler GUI/MultSegGLMEGUI.mlapp`
- `Analysis/Matlab Version/Tyler GUI/RWAnalysis2.m`
- `RW Task Patient Data/RWAnalysis2_MedianNorm_2sec.mat`

## Using GUI to Plot Fixed Effects

1. **Create an instance** of the `RWAnalysis2` class.  
2. **Set the path** to your `.mat` Analysis File:  
   `RWAnalysis2_MedianNorm_2sec.mat`  
3. **[Load](https://github.com/INMANLab/CAPTURE/blob/main/Analysis/Matlab%20Version/Tyler%20GUI/Docs/RWAnalysis2_Functions.md#loadmultdataobjvarargin) the multi-patient data:**  
   `GUIhandle.loadMultData;`  
4. **Launch the spectrogram GUI:**  
   `MultSegGLMEGUI(GUIhandle);`

Example code:
```matlab
GUIhandle = RWAnalysis2;
GUIhandle.AnalysisFile = '~/MATLAB/Data/RW Task Patient Data/RWAnalysis2_MedianNorm_2sec.mat';
GUIhandle.loadMultData;
MultSegGLMEGUI(GUIhandle);
```

5. **Input the desired parameters in the GUI:**

![GUI interface](https://github.com/INMANLab/CAPTURE/blob/main/Analysis/Matlab%20Version/Tyler%20GUI/Docs/MultSegGLMEGUI_interface.png)

### Core inputs
-  **Save Dir / Save Button**  
  Saves all generated figures and raw data as a `.mat` file with variables `MT` and `MTd` to the specified directory. Automatically creates the folder if it doesn't exist.

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

### GLME Model inputs

- **Dependent variable**  
  Selects the dependent variable of the model.  
  - ➤ This is the dropdown menu to the left of the '~' sign.  
  - ➤ Your selection will populate the 'GLME Model' text box.
  - ➤ For a description of each variable, refer to [getMultData function](https://github.com/INMANLab/CAPTURE/blob/main/Analysis/Matlab%20Version/Tyler%20GUI/Docs/RWAnalysis2_Functions.md#getmultdataobj-varargin)

- **Fixed or Random Effect selector**  
  Selects the Fixed Effect (FE) or Random Effect (RE) to be added to the model.  
  - ➤ This is the dropdown menu to the right of the '~' sign.
  - ➤ Selecting a variable and then pressing the `Add FE` button will add it as a Fixed Effect in the model. Pressing the `Add RE` will add it as a Random Effect.
 
- **Intercept checkbox**  
  Determines whether an intercept term is included in the GLME model.

  - ➤ If **checked**: The model includes a non-zero intercept term, allowing for a baseline level of the dependent variable when all predictors are zero.
  
  - ➤ If **unchecked**: The model **excludes** the intercept by prepending `-1` to the formula (e.g., `mtAlpha ~ -1 + Vel + Fix + ...`). This forces the model to go through the origin, assuming the baseline level of the dependent variable is zero when all predictors are zero.

- **GLME Model text box**  
  Can be populated manually or via the dropdown menu inputs above. Represents the string input for MATLAB's `fitglme` function.  
  - ➤ _NOTE:_ Using the dropdown menu inputs overrides any text manually inputted into the model. This is an issue that we need to fix.

  
The bottom of the GUI has the complete written out model that is trained on the data and serves as an example input. 
