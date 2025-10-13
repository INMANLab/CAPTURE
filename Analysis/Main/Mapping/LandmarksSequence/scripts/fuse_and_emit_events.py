#!/usr/bin/env python3
"""
Fuse beam and HMM results to emit final landmark events.

Strategy:
1. Load results from both beam filter and HMM
2. For each landmark, select the detection with higher confidence
3. Apply final quality gates (persistence, score threshold)
4. Emit final events.csv

Input:
    - out/events_beam.csv
    - out/events_hmm.csv
    - conf/sequence_priors.yaml

Output:
    - out/events.csv
"""
import sys
from pathlib import Path
import pandas as pd
import yaml
from typing import Dict


def load_sequence_config(config_path: Path) -> Dict:
    """Load sequence prior configuration."""
    with open(config_path, 'r') as f:
        return yaml.safe_load(f)


def fuse_events(beam_events: pd.DataFrame, hmm_events: pd.DataFrame, config: Dict) -> pd.DataFrame:
    """
    Fuse beam and HMM events to create final landmark events.

    Strategy:
    - For each landmark, prefer the detection with higher score
    - If landmark appears in only one method, use that detection
    - Apply quality gates from config

    Parameters
    ----------
    beam_events : pd.DataFrame
        Events from beam filter
    hmm_events : pd.DataFrame
        Events from HMM
    config : Dict
        Sequence configuration

    Returns
    -------
    pd.DataFrame
        Final fused events
    """
    print(f"\n{'='*70}")
    print("Fusing Beam and HMM Results")
    print(f"{'='*70}")

    # Get all unique landmarks from both sources
    all_landmarks = set(beam_events['landmark_id'].unique()) | set(hmm_events['landmark_id'].unique())
    all_landmarks = sorted(list(all_landmarks))

    fused_events = []

    for lm_id in all_landmarks:
        beam_det = beam_events[beam_events['landmark_id'] == lm_id]
        hmm_det = hmm_events[hmm_events['landmark_id'] == lm_id]

        if len(beam_det) == 0 and len(hmm_det) == 0:
            continue

        # Select best detection
        if len(beam_det) > 0 and len(hmm_det) > 0:
            # Both methods found this landmark - pick higher score
            beam_score = beam_det.iloc[0]['score']
            hmm_score = hmm_det.iloc[0]['score']

            if beam_score >= hmm_score:
                selected = beam_det.iloc[0]
                method = 'beam_preferred'
            else:
                selected = hmm_det.iloc[0]
                method = 'hmm_preferred'

            print(f"{lm_id}: Both methods - selected {method} (beam={beam_score:.3f}, hmm={hmm_score:.3f})")

        elif len(beam_det) > 0:
            selected = beam_det.iloc[0]
            method = 'beam_only'
            print(f"{lm_id}: Beam only (score={selected['score']:.3f})")

        else:
            selected = hmm_det.iloc[0]
            method = 'hmm_only'
            print(f"{lm_id}: HMM only (score={selected['score']:.3f})")

        fused_events.append({
            'landmark_id': selected['landmark_id'],
            'timestamp_sec': selected['timestamp_sec'],
            's_m': selected['s_m'],
            'score': selected['score'],
            'video_frame': selected['video_frame'],
            'walk': selected['walk'],
            'method': method
        })

    print(f"\n{'='*70}")
    print(f"Fused {len(fused_events)} landmark events")
    print(f"{'='*70}\n")

    return pd.DataFrame(fused_events)


def apply_quality_gates(events: pd.DataFrame, config: Dict) -> pd.DataFrame:
    """
    Apply final quality gates to events.

    Parameters
    ----------
    events : pd.DataFrame
        Fused events
    config : Dict
        Sequence configuration

    Returns
    -------
    pd.DataFrame
        Filtered events
    """
    min_score = config['gates']['min_score']

    print(f"\nApplying quality gates:")
    print(f"  Minimum score: {min_score}")

    # Filter by score
    n_before = len(events)
    events = events[events['score'] >= min_score]
    n_after = len(events)

    print(f"  Filtered: {n_before} → {n_after} events")

    return events


def main():
    """Main execution."""
    print("=" * 70)
    print("Fuse and Emit Events")
    print("=" * 70)

    # Load configurations
    seq_config = load_sequence_config(Path('conf/sequence_priors.yaml'))

    # Load beam and HMM results
    beam_path = Path('out/events_beam.csv')
    hmm_path = Path('out/events_hmm.csv')

    if not beam_path.exists():
        print(f"Error: Beam events not found at {beam_path}")
        print("Run ordered_beam_filter.py first")
        return

    if not hmm_path.exists():
        print(f"Error: HMM events not found at {hmm_path}")
        print("Run hmm_sequence_align.py first")
        return

    beam_events = pd.read_csv(beam_path)
    hmm_events = pd.read_csv(hmm_path)

    print(f"\nLoaded beam events: {len(beam_events)}")
    print(f"Loaded HMM events: {len(hmm_events)}")

    # Fuse events
    events = fuse_events(beam_events, hmm_events, seq_config)

    # Apply quality gates
    events = apply_quality_gates(events, seq_config)

    # Sort by timestamp
    events = events.sort_values('timestamp_sec').reset_index(drop=True)

    # Save final events
    output_path = Path('out/events.csv')
    events.to_csv(output_path, index=False)

    print(f"\n✓ Saved final events to: {output_path}")
    print(f"  Total events: {len(events)}")

    # Print summary
    print(f"\nFinal Event Summary:")
    print(f"  Total landmarks detected: {len(events)}")
    for _, event in events.iterrows():
        print(f"    {event['landmark_id']}: t={event['timestamp_sec']:.2f}s, "
              f"s={event['s_m']:.1f}m, score={event['score']:.3f} ({event['method']})")


if __name__ == "__main__":
    main()
