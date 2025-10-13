#!/usr/bin/env python3
"""
Build arc-length functions from XSens COM trajectory for multiple walks.

For each walk, extracts raw XSens CoM position from the Excel workbook
and computes cumulative distance (arc-length) along the trajectory.

This provides a dead-reckoned distance metric without requiring GPS alignment.

Usage:
    python build_arc_length_multi_walk.py --walks Walk1 Walk2 Walk3
    python build_arc_length_multi_walk.py --walks Walk1  # Single walk
    python build_arc_length_multi_walk.py --all          # Process all available walks

Input:
    - ~/Downloads/RW/RW1/{walk_name}/Xsens_*/RW_1_w*.xlsx

Output:
    - out/arc_length_{walk_name_lower}.csv (timestamp, s_m)
"""
import sys
from pathlib import Path
import pandas as pd
import numpy as np
import argparse

# Add parent to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent.parent))


def find_xsens_workbook(walk_name: str, base_dir: Path = None) -> Path:
    """
    Find XSens Excel workbook for a given walk.

    Parameters
    ----------
    walk_name : str
        Walk name (e.g., 'Walk1', 'Walk2')
    base_dir : Path, optional
        Base directory containing walks (default: ~/Downloads/RW/RW1)

    Returns
    -------
    Path
        Path to XSens workbook
    """
    if base_dir is None:
        base_dir = Path.home() / "Downloads/RW/RW1"

    walk_dir = base_dir / walk_name

    if not walk_dir.exists():
        raise FileNotFoundError(f"Walk directory not found: {walk_dir}")

    # Find Xsens folder
    xsens_folders = list(walk_dir.glob("Xsens*"))
    if not xsens_folders:
        raise FileNotFoundError(f"No Xsens folder found in {walk_dir}")

    xsens_folder = xsens_folders[0]

    # Find Excel workbook
    workbooks = list(xsens_folder.glob("*.xlsx"))
    if not workbooks:
        raise FileNotFoundError(f"No XSens workbook found in {xsens_folder}")

    return workbooks[0]


def extract_com_positions(workbook_path: Path) -> pd.DataFrame:
    """
    Extract Center of Mass positions from XSens workbook.

    Parameters
    ----------
    workbook_path : Path
        Path to XSens Excel workbook

    Returns
    -------
    pd.DataFrame
        DataFrame with columns: frame, com_x, com_y, com_z
    """
    print(f"  Reading workbook: {workbook_path.name}")

    # Read the Center of Mass sheet
    # XSens workbooks typically have a 'Center of Mass' sheet
    sheet_names = pd.ExcelFile(workbook_path).sheet_names
    print(f"  Available sheets: {sheet_names}")

    # Try to find COM sheet
    com_sheet = None
    for sheet in sheet_names:
        if 'center of mass' in sheet.lower() or 'com' in sheet.lower():
            com_sheet = sheet
            break

    if com_sheet is None:
        raise ValueError(f"Could not find Center of Mass sheet in {workbook_path}")

    print(f"  Reading sheet: {com_sheet}")
    df = pd.read_excel(workbook_path, sheet_name=com_sheet)

    print(f"  Sheet shape: {df.shape}")
    print(f"  Columns: {list(df.columns)[:10]}...")  # Show first 10 columns

    # Extract position columns
    # Common column names: 'Frame', 'CoM pos x', 'CoM pos y', 'CoM pos z'
    frame_col = None
    x_col = None
    y_col = None
    z_col = None

    for col in df.columns:
        col_lower = str(col).lower()
        if 'frame' in col_lower and frame_col is None:
            frame_col = col
        elif 'com pos x' in col_lower or (col_lower == 'x' and x_col is None):
            x_col = col
        elif 'com pos y' in col_lower or (col_lower == 'y' and y_col is None):
            y_col = col
        elif 'com pos z' in col_lower or (col_lower == 'z' and z_col is None):
            z_col = col

    if not all([frame_col, x_col, y_col, z_col]):
        raise ValueError(f"Could not find all required columns. Found: frame={frame_col}, x={x_col}, y={y_col}, z={z_col}")

    print(f"  Using columns: frame={frame_col}, x={x_col}, y={y_col}, z={z_col}")

    # Extract data
    com_df = df[[frame_col, x_col, y_col, z_col]].copy()
    com_df.columns = ['frame', 'com_x', 'com_y', 'com_z']

    # Remove NaN rows
    com_df = com_df.dropna()

    print(f"  Extracted {len(com_df)} COM positions")

    return com_df


def compute_arc_length(com_df: pd.DataFrame, walk_name: str, output_csv: Path):
    """
    Compute cumulative arc-length from CoM trajectory.

    Parameters
    ----------
    com_df : pd.DataFrame
        DataFrame with com_x, com_y, com_z columns
    walk_name : str
        Walk name for display
    output_csv : Path
        Path to output arc_length.csv
    """
    print(f"\n{'='*70}")
    print(f"Building Arc-Length for {walk_name}")
    print(f"{'='*70}")

    # Compute delta distance between consecutive points
    x = com_df['com_x'].values
    y = com_df['com_y'].values
    z = com_df['com_z'].values

    # Compute 3D Euclidean distance between consecutive samples
    dx = np.diff(x)
    dy = np.diff(y)
    dz = np.diff(z)
    d_dist = np.sqrt(dx**2 + dy**2 + dz**2)

    # Cumulative sum for arc-length
    # First sample is at s=0
    s_m = np.concatenate([[0.0], np.cumsum(d_dist)])

    # Create synthetic timestamp based on frame number
    # Assume ~100Hz sampling (typical for XSens)
    frame = com_df['frame'].values
    timestamp = frame / 100.0  # Approximate timestamp in seconds

    # Create output dataframe
    out_df = pd.DataFrame({
        'timestamp': timestamp,
        's_m': s_m
    })

    # Save
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(output_csv, index=False)

    print(f"\n✓ Saved arc-length to: {output_csv}")
    print(f"  Total samples: {len(out_df)}")
    print(f"  Arc-length range: {s_m.min():.2f} - {s_m.max():.2f} m")
    print(f"  Total distance: {s_m[-1]:.2f} m")
    print(f"  Time range: {timestamp.min():.2f} - {timestamp.max():.2f} s")
    print(f"  Duration: {(timestamp.max() - timestamp.min()) / 60:.1f} min")
    print(f"{'='*70}\n")


def process_walk(walk_name: str, output_dir: Path, base_dir: Path = None):
    """
    Process a single walk: extract COM and compute arc-length.

    Parameters
    ----------
    walk_name : str
        Walk name (e.g., 'Walk1', 'Walk2')
    output_dir : Path
        Output directory for arc_length CSVs
    base_dir : Path, optional
        Base directory containing walks

    Returns
    -------
    bool
        True if successful, False otherwise
    """
    print(f"\n{'='*70}")
    print(f"Processing {walk_name}")
    print(f"{'='*70}")

    try:
        # Find XSens workbook
        workbook_path = find_xsens_workbook(walk_name, base_dir)
        print(f"Found workbook: {workbook_path}")

        # Extract COM positions
        com_df = extract_com_positions(workbook_path)

        # Compute arc-length
        output_csv = output_dir / f"arc_length_{walk_name.lower()}.csv"
        compute_arc_length(com_df, walk_name, output_csv)

        return True

    except Exception as e:
        print(f"\n✗ Failed to process {walk_name}: {e}")
        import traceback
        traceback.print_exc()
        return False


def discover_walks(base_dir: Path = None):
    """
    Discover all available walks in the base directory.

    Parameters
    ----------
    base_dir : Path, optional
        Base directory to search (default: ~/Downloads/RW/RW1)

    Returns
    -------
    list
        List of walk names found
    """
    if base_dir is None:
        base_dir = Path.home() / "Downloads/RW/RW1"

    if not base_dir.exists():
        return []

    # Find all Walk* directories
    walk_dirs = sorted(base_dir.glob("Walk*"))
    return [d.name for d in walk_dirs if d.is_dir()]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Build arc-length functions from XSens data')
    parser.add_argument('--walks', nargs='+', help='Walk names to process (e.g., Walk1 Walk2)')
    parser.add_argument('--all', action='store_true', help='Process all available walks')
    parser.add_argument('--base-dir', type=str, help='Base directory containing walks (default: ~/Downloads/RW/RW1)')
    parser.add_argument('--output-dir', type=str, default='out', help='Output directory (default: out)')
    args = parser.parse_args()

    print("=" * 70)
    print("Building Arc-Length Functions for Multiple Walks")
    print("=" * 70)

    # Determine base directory
    base_dir = Path(args.base_dir).expanduser() if args.base_dir else None

    # Determine which walks to process
    if args.all:
        walks = discover_walks(base_dir)
        if not walks:
            print("No walks found!")
            sys.exit(1)
        print(f"\nDiscovered walks: {walks}")
    elif args.walks:
        walks = args.walks
    else:
        print("Error: Must specify --walks or --all")
        parser.print_help()
        sys.exit(1)

    output_dir = Path(args.output_dir)

    # Process all walks
    results = {}

    for walk_name in walks:
        success = process_walk(walk_name, output_dir, base_dir)
        results[walk_name] = success

    # Summary
    print("\n" + "=" * 70)
    print("Summary")
    print("=" * 70)
    for walk_name, success in results.items():
        status = "✓" if success else "✗"
        print(f"  {walk_name}: {status}")
    print("=" * 70)

    # Exit with error if any failed
    if not all(results.values()):
        sys.exit(1)
