# Python Pipeline	
# NeuroIOT Data: Tutorial & Usage Guide

### BEFORE DOING ANYTHING DOWNLOAD RW.zip FROM **SyncIOTDataPipeline** IN THE BOX -- Hopefully this moves to a server based file but for now everything must be downloaded locally.

## START HERE
Welcome to the NeuroIOT data repository! This tutorial will walk you through the process of extracting, merging, and using the various sensor data collected by Xsens, pupilphone(s), chestphone(s), and optional events.

Please read any sections marked **TUTORIAL** to set up your environment and understand how the provided scripts work. If you have any questions, refer to the **Troubleshooting** section below.

---

## TUTORIAL: Overview

1. **Xsens Data**  
   - Captured in `.xlsx` files with multiple sheets, each representing different sensor readings (e.g., orientation, acceleration, joint angles, etc.).
   - Use the script to specify which sheets to merge by modifying the `xsens_sheets` list.

2. **Pupilphone & Chestphone Data**  
   - Logged in CSV files, typically containing `[accelData, ambientLightData, gpsData, gyroData, magnetData]`.
   - You can choose which data types (e.g., `"accelData"`, `"gyroData"`) to include by editing the `pupilphone_sheets` and `chestphone_sheets` lists in the script.

3. **Events CSV**  
   - An additional CSV (e.g., `evnts_RWNApp_RW1_Walk1.csv`) that contains event timestamps under a column named **`NTP`**.
   - Set `include_events = True` (or `False`) to toggle whether to merge this file.

4. **Merging Logic**  
   - The script merges everything by matching on a master **`NTP_Timestamp`** column.
   - Rows that have an `NTP_Timestamp` but no other sensor data are automatically removed at the end.

---

## TUTORIAL: How to Use the Script

1. **Check Requirements**
   - Python 3.8+ (or newer).
   - Packages: `pandas`, `numpy` (implicit), `scipy` (for `.mat` files), and `datetime`.
   - Your `.mat` file must contain the correct NTP arrays (e.g., `ntp_xs`) that you want to use as a reference.

2. **Folder Structure**
   - A typical layout (based on the script’s defaults) might look like:

     ```
     C:\Users\burke\OneDrive\Desktop\
     └── RW
         └── RW1
             └── Walk1
                 ├── Xsens_1
                 │   └── (One or more Xsens .xlsx files)
                 ├── pupilphone_1
                 │   └── accelData.csv
                 │   └── gyroData.csv
                 ├── chestphone_1
                 │   └── accelData.csv
                 │   └── gyroData.csv
                 ├── RWNApp_RW1_Walk1.mat
                 └── evnts_RWNApp_RW1_Walk1.csv
     ```

   - Adjust file paths in the script as needed.

3. **Script Configuration**
   - Open the script in a text editor/IDE and locate the **Configuration** section:
     
     ```python
     # -----------------------------------------------------------------------------------
     # 📂 Configuration
     # -----------------------------------------------------------------------------------
     RW = 1
     walk = 1
     local = "C:\\Users\\burke\\OneDrive\\Desktop\\"

     ntp_source = "ntp_xs"  # or ntp_gp, ntp_np, etc.

     xsens_sheets = ["Center of Mass"]
     pupilphone_sheets = ["accelData", "gyroData"]
     chestphone_sheets = ["accelData", "gyroData"]

     include_ntp_pupil = False
     include_ntp_gp = False
     include_events = True
     ```
   - **`xsens_sheets`**: a list of strings matching the sheet names in the Xsens `.xlsx` file.  
   - **`pupilphone_sheets`** & **`chestphone_sheets`**: lists of data prefixes for CSV files in each phone’s folder.  
   - **`include_ntp_pupil`** & **`include_ntp_gp`**: toggles for whether to load pupil or GoPro NTP data from the `.mat` file.  
   - **`include_events`**: a boolean indicating whether to merge the events CSV.

4. **Run the Script**
   - Use the command line or an IDE (e.g., VSCode, PyCharm) to run your Python script:
     ```bash
     python merge_neuroiot.py
     ```
   - If everything is configured correctly, you should see messages about which files were found and a success message upon saving the final merged CSV:
     ```
     ✅ Merged Data Successfully Saved: C:\Users\burke\OneDrive\Desktop\RW\RW1\Walk1\Merged_Data_NTP.csv
     ```

5. **Output: `Merged_Data_NTP.csv`**
   - This file contains all timestamps in a single **`NTP_Timestamp`** column plus any data you selected from Xsens, pupilphone, chestphone, and events.
   - Any rows that contained **only** a timestamp and no other sensor/event data are automatically removed.

---

## 📂 Xsens Sheet Names (Available Options)

Below is a list of possible sheet names in the Xsens `.xlsx` file. Choose whichever you need by adding them to `xsens_sheets` in the script:

- `General Information`
- `Markers`
- `Segment Orientation - Quat`
- `Segment Orientation - Euler`
- `Segment Position`
- `Segment Velocity`
- `Segment Acceleration`
- `Segment Angular Velocity`
- `Segment Angular Acceleration`
- `Joint Angles ZXY`
- `Joint Angles XZY`
- `Ergonomic Joint Angles ZXY`
- `Ergonomic Joint Angles XZY`
- `Center of Mass`
- `Sensor Free Acceleration`
- `Sensor Magnetic Field`
- `Sensor Orientation - Quat`
- `Sensor Orientation - Euler`

---

## Pupils & Chestphone Data (Available CSV Prefixes)

For both pupilphone and chestphone, the following CSVs might exist:
- `accelData`
- `ambientLightData`
- `gpsData`
- `gyroData`
- `magnetData`

In the script, you specify these in `pupilphone_sheets` or `chestphone_sheets`. For instance:
```python
pupilphone_sheets = ["accelData", "gyroData"]
chestphone_sheets = ["accelData", "gyroData"]

