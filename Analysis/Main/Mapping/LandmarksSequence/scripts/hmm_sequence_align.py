#!/usr/bin/env python3
"""
HMM-based sequence alignment for landmark detection.

Uses a Hidden Markov Model with Viterbi decoding to find the most likely
sequence of landmark observations given:
1. Expected landmark order (route prior)
2. Arc-length distance priors
3. Vision scores
4. Transition probabilities (self-loop, forward, skip)

Input:
    - conf/route_prior.yaml
    - conf/sequence_priors.yaml
    - Analysis/Main/Mapping/Landmarks/out/detections/*.csv
    - out/arc_length.csv

Output:
    - out/events_hmm.csv
"""
import sys
from pathlib import Path
import pandas as pd
import numpy as np
import yaml
from typing import List, Dict, Tuple
from scipy.stats import norm

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

    Simplified approach: assumes temporal alignment between detection timestamps
    and arc-length samples.
    """
    detections = detections.copy()
    detections['s_m'] = np.nan

    for idx, row in detections.iterrows():
        timestamp_sec = row['timestamp_sec']
        arc_idx = int(timestamp_sec * 100)  # Assuming ~100Hz sampling
        if arc_idx < len(arc_length):
            detections.at[idx, 's_m'] = arc_length.iloc[arc_idx]['s_m']

    return detections


def build_emission_probs(detections: pd.DataFrame, landmarks: List[Dict], config: Dict) -> Tuple[np.ndarray, List]:
    """
    Build emission probability matrix P(observation | state).

    For each detection (observation), compute probability of being emitted by each
    landmark state based on:
    - Vision score (CLIP similarity)
    - Arc-length distance from expected position

    Parameters
    ----------
    detections : pd.DataFrame
        All detections with s_m and score columns
    landmarks : List[Dict]
        Ordered list of landmark info from route_prior
    config : Dict
        Sequence configuration

    Returns
    -------
    emission_probs : np.ndarray
        Shape (n_observations, n_states)
    observations : List
        List of detection indices
    """
    n_obs = len(detections)
    n_states = len(landmarks)

    emission_probs = np.zeros((n_obs, n_states))
    emission_floor = config['hmm']['emission_floor']

    weights = config['weights']
    w_vision = weights['w_vision']
    w_sdist = weights['w_sdist']

    for obs_idx, (det_idx, det) in enumerate(detections.iterrows()):
        det_lm_id = det['landmark_id']
        det_s = det['s_m']
        det_score = det['score']

        for state_idx, lm_info in enumerate(landmarks):
            lm_id = lm_info['id']
            s_expected = lm_info['s_m']
            s_std = lm_info['s_std']

            # Identity match: only consider detections for this landmark
            if det_lm_id != lm_id:
                emission_probs[obs_idx, state_idx] = emission_floor
                continue

            # Vision score probability (normalized)
            p_vision = det_score

            # Arc-length distance probability (Gaussian)
            p_sdist = norm.pdf(det_s, loc=s_expected, scale=s_std)
            # Normalize to [0, 1] range
            p_sdist = p_sdist / norm.pdf(s_expected, loc=s_expected, scale=s_std)

            # Combined probability
            p_emission = w_vision * p_vision + w_sdist * p_sdist

            emission_probs[obs_idx, state_idx] = max(p_emission, emission_floor)

    return emission_probs, list(detections.index)


def build_transition_matrix(n_states: int, config: Dict) -> np.ndarray:
    """
    Build state transition matrix A[i,j] = P(state_j | state_i).

    Transition rules:
    - Self-loop: stay in current landmark state
    - Forward: move to next landmark in sequence
    - Skip: skip one landmark and move to next+1

    Parameters
    ----------
    n_states : int
        Number of landmark states
    config : Dict
        Sequence configuration with HMM parameters

    Returns
    -------
    trans_matrix : np.ndarray
        Shape (n_states, n_states)
    """
    hmm_config = config['hmm']
    p_self = hmm_config['self_trans']
    p_fwd = hmm_config['fwd_trans']
    p_skip = hmm_config['skip_trans']

    trans_matrix = np.zeros((n_states, n_states))

    for i in range(n_states):
        # Self-loop
        trans_matrix[i, i] = p_self

        # Forward transition
        if i + 1 < n_states:
            trans_matrix[i, i + 1] = p_fwd

        # Skip transition
        if i + 2 < n_states:
            trans_matrix[i, i + 2] = p_skip

    # Normalize rows
    row_sums = trans_matrix.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1  # Avoid division by zero
    trans_matrix = trans_matrix / row_sums

    return trans_matrix


def viterbi(emission_probs: np.ndarray, trans_matrix: np.ndarray,
            start_probs: np.ndarray) -> Tuple[List[int], float]:
    """
    Viterbi algorithm for finding most likely state sequence.

    Parameters
    ----------
    emission_probs : np.ndarray
        Emission probabilities P(obs | state), shape (n_obs, n_states)
    trans_matrix : np.ndarray
        Transition probabilities P(state_j | state_i), shape (n_states, n_states)
    start_probs : np.ndarray
        Initial state probabilities, shape (n_states,)

    Returns
    -------
    path : List[int]
        Most likely state sequence (state indices)
    prob : float
        Log probability of best path
    """
    n_obs, n_states = emission_probs.shape

    # Viterbi tables
    viterbi_table = np.zeros((n_obs, n_states))
    backpointer = np.zeros((n_obs, n_states), dtype=int)

    # Initialize
    viterbi_table[0, :] = np.log(start_probs + 1e-10) + np.log(emission_probs[0, :] + 1e-10)

    # Forward pass
    for t in range(1, n_obs):
        for s in range(n_states):
            # Compute max probability of transitioning to state s at time t
            trans_probs = viterbi_table[t - 1, :] + np.log(trans_matrix[:, s] + 1e-10)
            max_prev_state = np.argmax(trans_probs)

            viterbi_table[t, s] = trans_probs[max_prev_state] + np.log(emission_probs[t, s] + 1e-10)
            backpointer[t, s] = max_prev_state

    # Backtrack
    path = [np.argmax(viterbi_table[-1, :])]
    for t in range(n_obs - 1, 0, -1):
        path.insert(0, backpointer[t, path[0]])

    best_prob = viterbi_table[-1, path[-1]]

    return path, best_prob


def hmm_align(detections: pd.DataFrame, route_prior: Dict, config: Dict) -> pd.DataFrame:
    """
    Apply HMM-based sequence alignment to select landmark events.

    Parameters
    ----------
    detections : pd.DataFrame
        All detections with s_m (arc-length) and score columns
    route_prior : Dict
        Route configuration with ordered landmarks
    config : Dict
        Sequence configuration

    Returns
    -------
    pd.DataFrame
        Selected landmark events
    """
    landmarks = route_prior['landmarks']
    n_states = len(landmarks)

    print(f"\n{'='*70}")
    print("HMM Sequence Alignment")
    print(f"{'='*70}")
    print(f"States (landmarks): {n_states}")
    print(f"Observations (detections): {len(detections)}")

    # Build emission probabilities
    emission_probs, obs_indices = build_emission_probs(detections, landmarks, config)

    # Build transition matrix
    trans_matrix = build_transition_matrix(n_states, config)

    # Initial state probabilities (uniform)
    start_probs = np.ones(n_states) / n_states

    # Run Viterbi
    print("\nRunning Viterbi decoding...")
    state_sequence, log_prob = viterbi(emission_probs, trans_matrix, start_probs)

    print(f"Best path log probability: {log_prob:.2f}")

    # Extract events: for each state, find best observation assigned to it
    events = []
    for state_idx, lm_info in enumerate(landmarks):
        lm_id = lm_info['id']

        # Find all observations assigned to this state
        obs_for_state = [obs_idx for obs_idx, s_idx in enumerate(state_sequence) if s_idx == state_idx]

        if not obs_for_state:
            print(f"\n{lm_id}: No observations assigned")
            continue

        # Select observation with highest emission probability
        best_obs_idx = max(obs_for_state, key=lambda o: emission_probs[o, state_idx])
        det_row = detections.iloc[best_obs_idx]

        print(f"\n{lm_id}: t={det_row['timestamp_sec']:.2f}s, "
              f"s={det_row['s_m']:.1f}m, score={det_row['score']:.3f}")

        events.append({
            'landmark_id': lm_id,
            'timestamp_sec': det_row['timestamp_sec'],
            's_m': det_row['s_m'],
            'score': det_row['score'],
            'video_frame': det_row['video_frame'],
            'walk': det_row['walk'],
            'method': 'hmm'
        })

    print(f"\n{'='*70}")
    print(f"HMM selected {len(events)} landmark events")
    print(f"{'='*70}\n")

    return pd.DataFrame(events)


def main():
    """Main execution."""
    print("=" * 70)
    print("HMM Sequence Alignment")
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

    # Apply minimum score filter
    min_score = seq_config['gates']['min_score']
    detections = detections[detections['score'] >= min_score]
    print(f"Detections above score threshold: {len(detections)}")

    # Sort by timestamp
    detections = detections.sort_values('timestamp_sec').reset_index(drop=True)

    # Apply HMM alignment
    events = hmm_align(detections, route_prior, seq_config)

    # Save results
    output_path = Path('out/events_hmm.csv')
    output_path.parent.mkdir(parents=True, exist_ok=True)
    events.to_csv(output_path, index=False)

    print(f"\n✓ Saved HMM events to: {output_path}")
    print(f"  Total events: {len(events)}")


if __name__ == "__main__":
    main()
