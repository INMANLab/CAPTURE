#!/usr/bin/env python3
"""
Ordered beam filter for landmark sequence detection.

Greedy baseline that uses:
1. Expected landmark order from route_prior.yaml
2. Arc-length distance windows
3. Vision score thresholding

Input:
    - conf/route_prior.yaml
    - conf/sequence_priors.yaml
    - Analysis/Main/Mapping/Landmarks/out/detections/*.csv (all detection files)
    - out/arc_length.csv

Output:
    - out/events_beam.csv
"""
import sys
from pathlib import Path
import pandas as pd
import numpy as np
import yaml
from typing import List, Dict

# Add parent to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent))


def load_route_prior(config_path: Path) -> Dict:
    """Load route prior configuration."""
    with open(config_path, 'r') as f:
        return yaml.safe_load(f)


def load_sequence_config(config_path: Path) -> Dict:
    """Load sequence prior configuration."""
    with open(config_path, 'r') as f:
        return yaml.safe_load(f)


def load_all_detections(detections_dir: Path) -> pd.DataFrame:
    """Load all detection CSV files and concatenate."""
    csv_files = list(detections_dir.glob("*.csv"))

    if not csv_files:
        raise ValueError(f"No detection CSV files found in {detections_dir}")

    dfs = []
    for csv_file in csv_files:
        df = pd.read_csv(csv_file)
        dfs.append(df)

    return pd.concat(dfs, ignore_index=True)


def merge_detections_with_arc_length(detections: pd.DataFrame, arc_length: pd.DataFrame) -> pd.DataFrame:
    """
    Merge detections with arc-length using nearest timestamp match.

    Parameters
    ----------
    detections : pd.DataFrame
        Detection data with timestamp_sec column
    arc_length : pd.DataFrame
        Arc-length data with timestamp and s_m columns

    Returns
    -------
    pd.DataFrame
        Merged dataframe with s_m column added
    """
    # Convert timestamp_sec to nanoseconds for matching
    # Assuming timestamp in arc_length is in nanoseconds and timestamp_sec is relative
    # We'll do a simple nearest neighbor join

    detections = detections.copy()
    detections['s_m'] = np.nan

    for idx, row in detections.iterrows():
        # Find nearest timestamp in arc_length
        # This is a simplified approach - in production you'd want proper time sync
        timestamp_sec = row['timestamp_sec']

        # For now, we'll use a simple approach:
        # Map timestamp_sec to arc_length by assuming linear progression
        # This is a placeholder - proper implementation would use actual time sync

        # Find the closest arc_length sample by index
        # Assuming detections are aligned with arc_length temporally
        arc_idx = int(timestamp_sec * 100)  # Assuming ~100Hz sampling
        if arc_idx < len(arc_length):
            detections.at[idx, 's_m'] = arc_length.iloc[arc_idx]['s_m']

    return detections


def beam_filter(detections: pd.DataFrame, route_prior: Dict, config: Dict) -> pd.DataFrame:
    """
    Apply ordered beam filter to select best landmark detections.

    Strategy:
    1. Process landmarks in order from route_prior
    2. For each landmark, find detections within arc-length window
    3. Select highest scoring detection above threshold
    4. Ensure temporal ordering (no going backwards in time)

    Parameters
    ----------
    detections : pd.DataFrame
        All detections with s_m (arc-length) column
    route_prior : Dict
        Route configuration with ordered landmarks
    config : Dict
        Sequence prior configuration

    Returns
    -------
    pd.DataFrame
        Filtered landmark events
    """
    landmarks = route_prior['landmarks']
    s_window = config['gates']['s_window_m']
    min_score = config['gates']['min_score']

    events = []
    last_timestamp = -np.inf

    print(f"\n{'='*70}")
    print("Beam Filter: Ordered Landmark Detection")
    print(f"{'='*70}")

    for lm_info in landmarks:
        lm_id = lm_info['id']
        s_expected = lm_info['s_m']
        s_std = lm_info['s_std']

        print(f"\n{lm_id}: Expected s={s_expected:.1f}m ± {s_std:.1f}m")

        # Filter detections for this landmark
        lm_dets = detections[detections['landmark_id'] == lm_id].copy()

        if len(lm_dets) == 0:
            print(f"  ⚠ No detections found")
            continue

        # Apply arc-length window
        s_min = s_expected - s_window
        s_max = s_expected + s_window
        lm_dets = lm_dets[(lm_dets['s_m'] >= s_min) & (lm_dets['s_m'] <= s_max)]

        if len(lm_dets) == 0:
            print(f"  ⚠ No detections in arc-length window [{s_min:.1f}, {s_max:.1f}]")
            continue

        # Filter by minimum score
        lm_dets = lm_dets[lm_dets['score'] >= min_score]

        if len(lm_dets) == 0:
            print(f"  ⚠ No detections above score threshold {min_score:.2f}")
            continue

        # Filter by temporal ordering (must be after last detection)
        lm_dets = lm_dets[lm_dets['timestamp_sec'] > last_timestamp]

        if len(lm_dets) == 0:
            print(f"  ⚠ No detections after last timestamp {last_timestamp:.2f}s")
            continue

        # Select highest scoring detection
        best_det = lm_dets.nlargest(1, 'score').iloc[0]

        print(f"  ✓ Selected: t={best_det['timestamp_sec']:.2f}s, "
              f"s={best_det['s_m']:.1f}m, score={best_det['score']:.3f}")

        events.append({
            'landmark_id': lm_id,
            'timestamp_sec': best_det['timestamp_sec'],
            's_m': best_det['s_m'],
            'score': best_det['score'],
            'video_frame': best_det['video_frame'],
            'walk': best_det['walk'],
            'method': 'beam'
        })

        last_timestamp = best_det['timestamp_sec']

    print(f"\n{'='*70}")
    print(f"Beam filter selected {len(events)} landmark events")
    print(f"{'='*70}\n")

    return pd.DataFrame(events)


def main():
    """Main execution."""
    print("=" * 70)
    print("Ordered Beam Filter")
    print("=" * 70)

    # Load configurations
    route_prior = load_route_prior(Path('conf/route_prior.yaml'))
    seq_config = load_sequence_config(Path('conf/sequence_priors.yaml'))

    # Load detections
    detections_dir = Path('../Landmarks/out')
    print(f"\nLoading detections from: {detections_dir}")
    detections = load_all_detections(detections_dir)
    print(f"Loaded {len(detections)} total detections")

    # Load arc-length
    arc_length_path = Path('out/arc_length.csv')
    print(f"Loading arc-length from: {arc_length_path}")
    arc_length = pd.read_csv(arc_length_path)
    print(f"Loaded {len(arc_length)} arc-length samples")

    # Merge detections with arc-length
    print("\nMerging detections with arc-length...")
    detections = merge_detections_with_arc_length(detections, arc_length)

    # Drop rows without arc-length match
    detections = detections.dropna(subset=['s_m'])
    print(f"Detections with arc-length: {len(detections)}")

    # Apply beam filter
    events = beam_filter(detections, route_prior, seq_config)

    # Save results
    output_path = Path('out/events_beam.csv')
    output_path.parent.mkdir(parents=True, exist_ok=True)
    events.to_csv(output_path, index=False)

    print(f"\n✓ Saved beam filter events to: {output_path}")
    print(f"  Total events: {len(events)}")


if __name__ == "__main__":
    main()
