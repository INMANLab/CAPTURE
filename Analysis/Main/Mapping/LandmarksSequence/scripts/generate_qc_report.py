#!/usr/bin/env python3
"""
Generate QC report for multi-walk landmark sequence recognition.

Evaluates:
- Sequence accuracy (detected vs expected landmarks)
- Arc-length alignment quality
- Score distribution
- Walk-to-walk consistency

Usage:
    python generate_qc_report.py
"""
import sys
from pathlib import Path
import pandas as pd
import numpy as np
import yaml
from typing import Dict

# Add parent to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent))


def load_route_prior(config_path: Path) -> Dict:
    """Load route prior configuration."""
    with open(config_path, 'r') as f:
        return yaml.safe_load(f)


def evaluate_sequence_accuracy(events: pd.DataFrame, route_prior: Dict) -> Dict:
    """
    Evaluate how well the detected sequence matches expected landmarks.
    """
    landmarks = route_prior['landmarks']
    lm_ids = [lm['id'] for lm in landmarks]

    results = {
        'total_expected': len(lm_ids),
        'total_detected': 0,
        'per_walk': {},
        'missing_landmarks': [],
        'detected_landmarks': []
    }

    for walk in events['walk'].unique():
        walk_events = events[events['walk'] == walk]
        detected_ids = set(walk_events['landmark_id'].unique())
        expected_ids = set(lm_ids)

        n_detected = len(detected_ids & expected_ids)
        n_expected = len(expected_ids)

        results['per_walk'][walk] = {
            'detected': n_detected,
            'expected': n_expected,
            'accuracy': n_detected / n_expected,
            'missing': list(expected_ids - detected_ids),
            'detected_list': list(detected_ids)
        }

    # Overall statistics
    all_detected = set(events['landmark_id'].unique())
    results['total_detected'] = len(all_detected)
    results['missing_landmarks'] = list(set(lm_ids) - all_detected)
    results['detected_landmarks'] = list(all_detected)

    return results


def evaluate_arc_length_alignment(events: pd.DataFrame, route_prior: Dict) -> Dict:
    """
    Evaluate how well detected arc-lengths match expected positions.
    """
    landmarks = {lm['id']: lm for lm in route_prior['landmarks']}

    results = {
        'per_landmark': {},
        'mean_abs_error': 0.0,
        'rmse': 0.0
    }

    errors = []

    for lm_id in events['landmark_id'].unique():
        lm_events = events[events['landmark_id'] == lm_id]
        expected_s = landmarks[lm_id]['s_m']

        detected_s = lm_events['s_m'].values
        errors_lm = detected_s - expected_s

        results['per_landmark'][lm_id] = {
            'expected_s': expected_s,
            'detected_s_mean': float(np.mean(detected_s)),
            'detected_s_std': float(np.std(detected_s)),
            'error_mean': float(np.mean(errors_lm)),
            'error_std': float(np.std(errors_lm)),
            'n_detections': len(detected_s)
        }

        errors.extend(errors_lm)

    if errors:
        results['mean_abs_error'] = float(np.mean(np.abs(errors)))
        results['rmse'] = float(np.sqrt(np.mean(np.square(errors))))

    return results


def evaluate_score_distribution(events: pd.DataFrame) -> Dict:
    """
    Evaluate vision score distribution.
    """
    scores = events['score'].values

    return {
        'mean': float(np.mean(scores)),
        'std': float(np.std(scores)),
        'min': float(np.min(scores)),
        'max': float(np.max(scores)),
        'median': float(np.median(scores)),
        'q25': float(np.percentile(scores, 25)),
        'q75': float(np.percentile(scores, 75))
    }


def generate_markdown_report(events: pd.DataFrame, route_prior: Dict, output_path: Path):
    """
    Generate markdown QC report.
    """
    # Evaluate metrics
    seq_accuracy = evaluate_sequence_accuracy(events, route_prior)
    arc_alignment = evaluate_arc_length_alignment(events, route_prior)
    score_dist = evaluate_score_distribution(events)

    # Generate report
    report_lines = []

    # Header
    report_lines.append("# Landmark Sequence Recognition QC Report")
    report_lines.append("")
    report_lines.append(f"**Generated:** {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}")
    report_lines.append("")

    # Executive Summary
    report_lines.append("## Executive Summary")
    report_lines.append("")
    report_lines.append(f"- **Total Events:** {len(events)}")
    report_lines.append(f"- **Walks Processed:** {events['walk'].nunique()}")
    report_lines.append(f"- **Landmarks Detected:** {seq_accuracy['total_detected']}/{seq_accuracy['total_expected']}")
    report_lines.append(f"- **Mean Arc-Length Error:** {arc_alignment['mean_abs_error']:.1f} m")
    report_lines.append(f"- **RMSE:** {arc_alignment['rmse']:.1f} m")
    report_lines.append(f"- **Mean Detection Score:** {score_dist['mean']:.3f}")
    report_lines.append("")

    # Per-Walk Analysis
    report_lines.append("## Per-Walk Analysis")
    report_lines.append("")

    for walk in sorted(events['walk'].unique()):
        walk_data = seq_accuracy['per_walk'][walk]
        report_lines.append(f"### {walk}")
        report_lines.append("")
        report_lines.append(f"- **Detected Landmarks:** {walk_data['detected']}/{walk_data['expected']} ({walk_data['accuracy']*100:.1f}%)")

        if walk_data['missing']:
            report_lines.append(f"- **Missing:** {', '.join(walk_data['missing'])}")

        walk_events = events[events['walk'] == walk]
        report_lines.append(f"- **Score Range:** [{walk_events['score'].min():.3f}, {walk_events['score'].max():.3f}]")
        report_lines.append("")

    # Landmark-Level Analysis
    report_lines.append("## Landmark-Level Analysis")
    report_lines.append("")
    report_lines.append("| Landmark | Expected s (m) | Detected s (m) | Error (m) | Detections | Score |")
    report_lines.append("|----------|----------------|----------------|-----------|------------|-------|")

    for lm_id in sorted(arc_alignment['per_landmark'].keys()):
        lm_data = arc_alignment['per_landmark'][lm_id]
        lm_events = events[events['landmark_id'] == lm_id]

        score_mean = lm_events['score'].mean()

        report_lines.append(
            f"| {lm_id} | "
            f"{lm_data['expected_s']:.1f} | "
            f"{lm_data['detected_s_mean']:.1f} ± {lm_data['detected_s_std']:.1f} | "
            f"{lm_data['error_mean']:.1f} ± {lm_data['error_std']:.1f} | "
            f"{lm_data['n_detections']} | "
            f"{score_mean:.3f} |"
        )

    report_lines.append("")

    # Missing Landmarks
    if seq_accuracy['missing_landmarks']:
        report_lines.append("## Missing Landmarks")
        report_lines.append("")
        report_lines.append("The following landmarks were not detected in any walk:")
        report_lines.append("")
        for lm_id in seq_accuracy['missing_landmarks']:
            lm_info = next(lm for lm in route_prior['landmarks'] if lm['id'] == lm_id)
            report_lines.append(f"- **{lm_id}**: {lm_info['name']} (expected at {lm_info['s_m']:.1f}m)")
        report_lines.append("")

    # Score Distribution
    report_lines.append("## Detection Score Distribution")
    report_lines.append("")
    report_lines.append(f"- **Mean:** {score_dist['mean']:.3f}")
    report_lines.append(f"- **Median:** {score_dist['median']:.3f}")
    report_lines.append(f"- **Std Dev:** {score_dist['std']:.3f}")
    report_lines.append(f"- **Range:** [{score_dist['min']:.3f}, {score_dist['max']:.3f}]")
    report_lines.append(f"- **IQR:** [{score_dist['q25']:.3f}, {score_dist['q75']:.3f}]")
    report_lines.append("")

    # QA Gates
    report_lines.append("## QA Gates")
    report_lines.append("")

    # Compute overall accuracy
    total_possible = seq_accuracy['total_expected'] * events['walk'].nunique()
    total_detected = len(events)
    overall_accuracy = total_detected / total_possible

    gates = [
        ("Sequence Coverage", overall_accuracy >= 0.80, f"{overall_accuracy:.2%} {'✓' if overall_accuracy >= 0.80 else '✗'}"),
        ("Mean Arc-Length Error", arc_alignment['mean_abs_error'] <= 100, f"{arc_alignment['mean_abs_error']:.1f}m {'✓' if arc_alignment['mean_abs_error'] <= 100 else '✗'}"),
        ("Mean Detection Score", score_dist['mean'] >= 0.28, f"{score_dist['mean']:.3f} {'✓' if score_dist['mean'] >= 0.28 else '✗'}")
    ]

    for gate_name, passed, value in gates:
        status = "✅ PASS" if passed else "❌ FAIL"
        report_lines.append(f"- **{gate_name}:** {status} ({value})")

    report_lines.append("")

    # Overall Assessment
    all_passed = all(g[1] for g in gates)
    report_lines.append("## Overall Assessment")
    report_lines.append("")
    if all_passed:
        report_lines.append("**✅ ALL QA GATES PASSED**")
        report_lines.append("")
        report_lines.append("The landmark sequence recognition pipeline meets quality standards.")
    else:
        report_lines.append("**⚠️ SOME QA GATES FAILED**")
        report_lines.append("")
        report_lines.append("Review failed gates and consider:")
        report_lines.append("- Adjusting detection thresholds")
        report_lines.append("- Refining exemplar images")
        report_lines.append("- Reviewing landmark definitions")

    report_lines.append("")

    # Write report
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, 'w') as f:
        f.write('\n'.join(report_lines))


def main():
    """Main execution."""
    print("=" * 70)
    print("Generating QC Report")
    print("=" * 70)

    # Load data
    route_prior = load_route_prior(Path('conf/route_prior.yaml'))
    events_path = Path('out/events_beam_multiwalk.csv')

    if not events_path.exists():
        print(f"\nError: Events file not found: {events_path}")
        print("Run sequence pipeline first!")
        sys.exit(1)

    events = pd.read_csv(events_path)
    print(f"\nLoaded {len(events)} events from {events['walk'].nunique()} walks")

    # Generate report
    output_path = Path('out/qc/sequence_qc_report.md')
    generate_markdown_report(events, route_prior, output_path)

    print(f"\n✓ QC report saved to: {output_path}")


if __name__ == "__main__":
    main()
