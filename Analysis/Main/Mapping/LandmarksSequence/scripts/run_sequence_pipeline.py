#!/usr/bin/env python3
"""
Run complete sequence recognition pipeline for multiple walks.

Processes each walk independently, then aggregates results.

Usage:
    python run_sequence_pipeline.py --walks Walk1 Walk2 Walk3 Walk4
    python run_sequence_pipeline.py --all
"""
import sys
from pathlib import Path
import pandas as pd
import numpy as np
import yaml
import argparse
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


def load_detection_for_walk(detections_dir: Path, walk_name: str) -> pd.DataFrame:
    """Load detection CSV for a specific walk."""
    # Find CSV file matching walk name
    pattern = f"{walk_name.lower()}_*.csv"
    csv_files = list(detections_dir.glob(pattern))

    if not csv_files:
        raise ValueError(f"No detection file found for {walk_name} in {detections_dir}")

    if len(csv_files) > 1:
        print(f"Warning: Multiple files found for {walk_name}, using {csv_files[0].name}")

    df = pd.read_csv(csv_files[0])

    # Ensure walk column exists
    if 'walk' not in df.columns:
        df['walk'] = walk_name

    return df


def merge_detections_with_arc_length(detections: pd.DataFrame, arc_length: pd.DataFrame, walk_name: str) -> pd.DataFrame:
    """
    Merge detections with arc-length using timestamp alignment.

    Uses approximate timestamp matching based on video FPS.
    """
    detections = detections.copy()
    detections['s_m'] = np.nan

    # For each detection, find nearest arc-length sample
    for idx, row in detections.iterrows():
        timestamp_sec = row['timestamp_sec']

        # Assuming ~100Hz XSens sampling
        arc_idx = int(timestamp_sec * 100)

        if arc_idx < len(arc_length):
            detections.at[idx, 's_m'] = arc_length.iloc[arc_idx]['s_m']

    return detections


def beam_filter(detections: pd.DataFrame, route_prior: Dict, config: Dict, walk_name: str) -> pd.DataFrame:
    """
    Apply ordered beam filter for single walk.
    """
    landmarks = route_prior['landmarks']
    s_window = config['gates']['s_window_m']
    min_score = config['gates']['min_score']

    events = []
    last_timestamp = -np.inf

    print(f"\n{'='*70}")
    print(f"Beam Filter: {walk_name}")
    print(f"{'='*70}")

    for lm_info in landmarks:
        lm_id = lm_info['id']
        s_expected = lm_info['s_m']
        s_std = lm_info['s_std']

        print(f"\n{lm_id}: Expected s={s_expected:.1f}m ± {s_std:.1f}m")

        # Filter detections for this landmark
        lm_dets = detections[detections['landmark_id'] == lm_id].copy()

        if len(lm_dets) == 0:
            print(f"  ⚠ No detections")
            continue

        # Apply arc-length window
        s_min = s_expected - s_window
        s_max = s_expected + s_window
        lm_dets = lm_dets[(lm_dets['s_m'] >= s_min) & (lm_dets['s_m'] <= s_max)]

        if len(lm_dets) == 0:
            print(f"  ⚠ No detections in window [{s_min:.1f}, {s_max:.1f}]")
            continue

        # Filter by score
        lm_dets = lm_dets[lm_dets['score'] >= min_score]

        if len(lm_dets) == 0:
            print(f"  ⚠ No detections above threshold {min_score:.2f}")
            continue

        # Temporal ordering
        lm_dets = lm_dets[lm_dets['timestamp_sec'] > last_timestamp]

        if len(lm_dets) == 0:
            print(f"  ⚠ No detections after t={last_timestamp:.2f}s")
            continue

        # Select best
        best_det = lm_dets.nlargest(1, 'score').iloc[0]

        print(f"  ✓ t={best_det['timestamp_sec']:.2f}s, "
              f"s={best_det['s_m']:.1f}m, score={best_det['score']:.3f}")

        events.append({
            'landmark_id': lm_id,
            'timestamp_sec': best_det['timestamp_sec'],
            's_m': best_det['s_m'],
            'score': best_det['score'],
            'video_frame': best_det['video_frame'],
            'walk': walk_name,
            'method': 'beam'
        })

        last_timestamp = best_det['timestamp_sec']

    print(f"\n{'='*70}")
    print(f"{walk_name}: {len(events)} events")
    print(f"{'='*70}\n")

    return pd.DataFrame(events)


def process_walk(walk_name: str, route_prior: Dict, seq_config: Dict,
                 detections_dir: Path, arc_length_dir: Path) -> pd.DataFrame:
    """
    Process single walk through sequence pipeline.
    """
    print(f"\n{'#'*70}")
    print(f"Processing {walk_name}")
    print(f"{'#'*70}")

    # Load detection data
    detections = load_detection_for_walk(detections_dir, walk_name)
    print(f"Loaded {len(detections)} detections for {walk_name}")

    # Load arc-length data
    arc_length_file = arc_length_dir / f"arc_length_{walk_name.lower()}.csv"
    if not arc_length_file.exists():
        print(f"Error: Arc-length file not found: {arc_length_file}")
        return pd.DataFrame()

    arc_length = pd.read_csv(arc_length_file)
    print(f"Loaded {len(arc_length)} arc-length samples")

    # Merge with arc-length
    detections = merge_detections_with_arc_length(detections, arc_length, walk_name)
    detections = detections.dropna(subset=['s_m'])
    print(f"Detections with arc-length: {len(detections)}")

    # Apply beam filter
    events = beam_filter(detections, route_prior, seq_config, walk_name)

    return events


def main():
    """Main execution."""
    parser = argparse.ArgumentParser(description='Run sequence pipeline on multiple walks')
    parser.add_argument('--walks', nargs='+', help='Walk names (e.g., Walk1 Walk2)')
    parser.add_argument('--all', action='store_true', help='Process all available walks')
    args = parser.parse_args()

    print("=" * 70)
    print("Multi-Walk Sequence Recognition Pipeline")
    print("=" * 70)

    # Load configurations
    route_prior = load_route_prior(Path('conf/route_prior.yaml'))
    seq_config = load_sequence_config(Path('conf/sequence_priors.yaml'))

    # Determine walks to process
    detections_dir = Path('../Landmarks/out/detections')
    arc_length_dir = Path('out')

    if args.all:
        # Discover all walks with both detections and arc-length
        detection_files = list(detections_dir.glob('walk*_*.csv'))
        walks = sorted(set([f.name.split('_')[0].capitalize() for f in detection_files]))
        print(f"\nDiscovered walks: {walks}")
    elif args.walks:
        walks = args.walks
    else:
        print("Error: Must specify --walks or --all")
        parser.print_help()
        sys.exit(1)

    # Process each walk
    all_events = []

    for walk_name in walks:
        try:
            events = process_walk(walk_name, route_prior, seq_config,
                                 detections_dir, arc_length_dir)
            if len(events) > 0:
                all_events.append(events)
        except Exception as e:
            print(f"\n✗ Error processing {walk_name}: {e}")
            import traceback
            traceback.print_exc()
            continue

    # Combine results
    if all_events:
        combined_events = pd.concat(all_events, ignore_index=True)

        # Save combined results
        output_path = Path('out/events_beam_multiwalk.csv')
        output_path.parent.mkdir(parents=True, exist_ok=True)
        combined_events.to_csv(output_path, index=False)

        print(f"\n{'='*70}")
        print("Pipeline Complete")
        print(f"{'='*70}")
        print(f"\n✓ Saved combined events to: {output_path}")
        print(f"  Total events: {len(combined_events)}")
        print(f"  Walks processed: {combined_events['walk'].nunique()}")

        # Summary per walk
        print("\nEvents per walk:")
        for walk in combined_events['walk'].unique():
            walk_events = combined_events[combined_events['walk'] == walk]
            print(f"  {walk}: {len(walk_events)} events")
    else:
        print("\n✗ No events generated")
        sys.exit(1)


if __name__ == "__main__":
    main()
