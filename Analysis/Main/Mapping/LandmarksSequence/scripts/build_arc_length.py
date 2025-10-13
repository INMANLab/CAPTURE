#!/usr/bin/env python3
"""
Build arc-length function from XSens COM trajectory.

Computes cumulative distance ŝ(t) along the route by integrating
dead-reckoned velocity from XSens data.

Input:
    - Analysis/Main/Mapping/data/xsens_com_enu.csv

Output:
    - out/arc_length.csv (timestamp, s_m)
"""
import sys
from pathlib import Path
import pandas as pd
import numpy as np

# Add parent to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent.parent))


def compute_arc_length(xsens_csv: Path, output_csv: Path):
    """
    Compute cumulative arc-length from XSens COM trajectory.

    Parameters
    ----------
    xsens_csv : Path
        Path to xsens_com_enu.csv
    output_csv : Path
        Path to output arc_length.csv
    """
    print("=" * 70)
    print("Building Arc-Length Function")
    print("=" * 70)
    print(f"\nInput: {xsens_csv}")

    # Load XSens data
    df = pd.read_csv(xsens_csv)
    print(f"Loaded {len(df)} XSens samples")
    print(f"Time range: {df['timestamp'].min():.2f} - {df['timestamp'].max():.2f}")

    # Compute delta distance between consecutive points
    # Using East-North position
    east = df['east_m'].values
    north = df['north_m'].values

    # Compute Euclidean distance between consecutive samples
    d_east = np.diff(east)
    d_north = np.diff(north)
    d_dist = np.sqrt(d_east**2 + d_north**2)

    # Cumulative sum for arc-length
    # First sample is at s=0
    s_m = np.concatenate([[0.0], np.cumsum(d_dist)])

    # Create output dataframe
    out_df = pd.DataFrame({
        'timestamp': df['timestamp'].values,
        's_m': s_m
    })

    # Save
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(output_csv, index=False)

    print(f"\n✓ Saved arc-length to: {output_csv}")
    print(f"  Total samples: {len(out_df)}")
    print(f"  Arc-length range: {s_m.min():.2f} - {s_m.max():.2f} m")
    print(f"  Total distance: {s_m[-1]:.2f} m")
    print("=" * 70)


if __name__ == "__main__":
    # Paths - use relative path from LandmarksSequence directory
    xsens_csv = Path("../data/xsens_com_enu.csv")
    output_csv = Path("out/arc_length.csv")

    if not xsens_csv.exists():
        print(f"Error: XSens data not found at {xsens_csv}")
        sys.exit(1)

    compute_arc_length(xsens_csv, output_csv)
