import os
import numpy as np
import pandas as pd
from datetime import datetime, timezone          # ← timezone added
import pytz                                       # ← pytz added
import scipy.io

pacific = pytz.timezone("America/Los_Angeles")

# -----------------------------------------------------------------------------------
# 📂 Configuration
# -----------------------------------------------------------------------------------
RW   = 1
walk = 2
local = "C:\\Users\\burke\\OneDrive\\Desktop\\"

ntp_source        = "ntp_xs"
neuralpace_enable = True

xsens_sheets       = []
pupilphone_sheets  = []
chestphone_sheets  = []
pupil_sheets       = ["gaze"]

include_ntp_pupil  = True
include_ntp_gp     = False
include_events     = True

fileName   = "Merged_example"
output_csv = f"{local}RW\\RW{RW}\\Walk{walk}\\{fileName}.csv"

# -----------------------------------------------------------------------------------
# ✅ NOTHING BELOW THIS LINE CHANGES ANY LOGIC
# -----------------------------------------------------------------------------------
xs_folder         = f"{local}RW\\RW{RW}\\Walk{walk}\\Xsens_{walk}\\"
mat_file_path     = f"{local}RW\\RW{RW}\\Walk{walk}\\RWNApp_RW{RW}_Walk{walk}.mat"
pupilphone_folder = f"{local}RW\\RW{RW}\\Walk{walk}\\pupilphone_{walk}\\"
chestphone_folder = f"{local}RW\\RW{RW}\\Walk{walk}\\chestphone_{walk}\\"
pupil_folder      = f"{local}RW\\RW{RW}\\Walk{walk}\\Pupil_{walk}\\"
events_file       = f"{local}RW\\RW{RW}\\Walk{walk}\\evnts_RWNApp_RW{RW}_Walk{walk}.csv"

# ---------------------------------------------------------------- check Xsens ------
xlsx_files = [f for f in os.listdir(xs_folder) if f.endswith(".xlsx")]
if not xlsx_files:
    print(f"❌ No .xlsx in {xs_folder}")
    raise SystemExit
xsens_file = os.path.join(xs_folder, xlsx_files[0])

# ---------------------------------------------------------------- load .mat --------
try:
    mat = scipy.io.loadmat(mat_file_path)
    if ntp_source not in mat:
        print(f"❌ '{ntp_source}' not found in .mat")
        raise SystemExit
    ntp_list = mat[ntp_source].flatten().tolist()

    ntp_pupil_data = (
        pd.DataFrame({
            "NTP_Pupil_Timestamp": mat["ntp_pupil"].flatten(),
            "Frame_Pupil": range(len(mat["ntp_pupil"]))
        }) if include_ntp_pupil and "ntp_pupil" in mat else pd.DataFrame()
    )

    ntp_gp_data = (
        pd.DataFrame({
            "NTP_GP_Timestamp": mat["ntp_gp"].flatten(),
            "Frame_GP": range(len(mat["ntp_gp"]))
        }) if include_ntp_gp and "ntp_gp" in mat else pd.DataFrame()
    )

    if neuralpace_enable:
        if "d_np" not in mat or "ntp_np" not in mat:
            print("❌ 'd_np' or 'ntp_np' missing in .mat")
            raise SystemExit
        d_np   = np.asarray(mat["d_np"])
        ntp_np = np.asarray(mat["ntp_np"]).flatten()
        neuralpace_cols = [f"NeuralPace_Ch{i+1}" for i in range(d_np.shape[1])]
        neuralpace_data = pd.DataFrame(d_np, columns=neuralpace_cols)
        neuralpace_data["NTP_Timestamp"] = ntp_np

except Exception as e:
    print(f"❌ .mat read error: {e}")
    raise SystemExit

# -------------------------------------------------------- load Xsens workbook ------
try:
    excel_data = pd.ExcelFile(xsens_file)
    xsens_data = pd.DataFrame({"NTP_Timestamp": ntp_list[:len(ntp_list)]})
    for sheet in xsens_sheets:
        if sheet not in excel_data.sheet_names:
            print(f"⚠ Sheet '{sheet}' missing – skipped.")
            continue
        df = excel_data.parse(sheet)
        df.columns = [f"{sheet}_{c}" for c in df.columns]
        xsens_data = pd.concat([xsens_data, df.iloc[:len(ntp_list)]], axis=1)
except Exception as e:
    print(f"❌ Xsens .xlsx error: {e}")
    raise SystemExit

# ------------------------------------------------------- timestamp converter -------
MATLAB_OFFSET_S = 62_167_219_200  # 0000-01-01 → 1970-01-01

def convert_to_matlab_ntp(ts):
    """
    • If ts > 1e12  → treat as nanoseconds and run the exact Pacific-time logic.
    • If ts is already a big MATLAB-seconds number → return unchanged.
    • Digit strings handled the same way.
    • Date-strings ('YYYY-MM-DD_HH-MM-SS-fff') fall back to original parser.
    """
    try:
        # ---------- numeric input -------------------------------------------------
        if isinstance(ts, (int, float, np.integer, np.floating)):
            val = float(ts)
            if val > 1e12:                       # nanoseconds
                timestamp_sec = val / 1e9
                dt_utc   = datetime.fromtimestamp(timestamp_sec, tz=timezone.utc)
                dt_local = dt_utc.astimezone(pacific)

                days_since_0000 = dt_local.toordinal() + 366
                frac = (dt_local.hour*3600 + dt_local.minute*60 +
                        dt_local.second + dt_local.microsecond/1e6) / 86400
                return (days_since_0000 + frac) * 86400

            return val                            # already big seconds

        # ---------- digit-string input -------------------------------------------
        s = str(ts).strip()
        if s.isdigit():
            return convert_to_matlab_ntp(float(s))

        # ---------- fallback: date-string -----------------------------------------
        matlab_epoch = datetime(1, 1, 1)
        dt = datetime.strptime(s, "%Y-%m-%d_%H-%M-%S-%f")
        dn = ((dt - matlab_epoch).total_seconds() + 172800) / 86400 + 365
        return dn * 86400

    except Exception:
        return None

# ------------------------------------------------------ generic CSV loader ---------
def load_sensor_data(folder, prefixes):
    if not prefixes:
        return pd.DataFrame(columns=["NTP_Timestamp"])
    files = [f for f in os.listdir(folder)
             if f.endswith(".csv") and any(f.startswith(p) for p in prefixes)]
    if not files:
        return pd.DataFrame(columns=["NTP_Timestamp"])

    out = pd.DataFrame()
    for f in files:
        path = os.path.join(folder, f)
        try:
            if "gaze" in f.lower():
                df = pd.read_csv(path, low_memory=False)
                df = df.iloc[:, 2:]                       # drop first 2 cols
                df.rename(columns={df.columns[0]: "Timestamp"}, inplace=True)
            else:
                df = pd.read_csv(path, header=None, dtype={0: str}, low_memory=False)
                base = os.path.splitext(f)[0]
                df.columns = ["Timestamp"] + [f"{base}_{i}" for i in range(1, len(df.columns))]

            df["NTP_Timestamp"] = df["Timestamp"].apply(convert_to_matlab_ntp)
            df.drop(columns=["Timestamp"], inplace=True)
            out = df if out.empty else pd.merge(out, df, on="NTP_Timestamp", how="outer")

        except Exception as e:
            print(f"❌ CSV error {f}: {e}")
    return out

# ------------------------------------------------------- load sensor folders -------
pupilphone_data = load_sensor_data(pupilphone_folder, pupilphone_sheets)
chestphone_data = load_sensor_data(chestphone_folder, chestphone_sheets)
pupil_data      = load_sensor_data(pupil_folder, pupil_sheets)

all_phone = pupil_data.merge(pupilphone_data, on="NTP_Timestamp", how="outer") \
                      .merge(chestphone_data, on="NTP_Timestamp", how="outer")

# -------------------------------------------------------------- merge everything --
merged = xsens_data.merge(all_phone, on="NTP_Timestamp", how="outer")
if neuralpace_enable:
    merged = merged.merge(neuralpace_data, on="NTP_Timestamp", how="outer")

if not ntp_pupil_data.empty:
    ntp_pupil_data.rename(columns={"NTP_Pupil_Timestamp": "NTP_Timestamp"}, inplace=True)
    merged = merged.merge(ntp_pupil_data, on="NTP_Timestamp", how="outer")
if not ntp_gp_data.empty:
    ntp_gp_data.rename(columns={"NTP_GP_Timestamp": "NTP_Timestamp"}, inplace=True)
    merged = merged.merge(ntp_gp_data, on="NTP_Timestamp", how="outer")

# --------------------------------------------------------------- events -----------
if include_events:
    try:
        ev = pd.read_csv(events_file)
        if "NTP" not in ev.columns:
            raise ValueError("'NTP' column missing in events CSV")
        ev.rename(columns={"NTP": "NTP_Timestamp"}, inplace=True)
        merged = merged.merge(ev, on="NTP_Timestamp", how="outer")
        print("✅ Events merged")
    except Exception as e:
        print(f"❌ Events merge error: {e}")
else:
    print("⚠ Events merge skipped")

# ------------------------------------------------ final clean & save --------------
merged = merged.dropna(axis=0, how="all", subset=merged.columns[1:])
merged.to_csv(output_csv, index=False)
print(f"✅ Merged data saved → {output_csv}")
