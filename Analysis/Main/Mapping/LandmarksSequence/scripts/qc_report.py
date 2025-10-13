#!/usr/bin/env python3
"""
Quality control report for landmark sequence detection.

Evaluates:
1. Sequence accuracy: Are landmarks in correct order?
2. Detection rate: How many landmarks were found?
3. Score distributions
4. Arc-length alignment quality

Input:
    - out/events.csv
    - conf/route_prior.yaml
    - (optional) ground truth file for evaluation

Output:
    - out/qc/sequence_qc_report.md
"""
import sys
from pathlib import Path
import pandas as pd
import numpy as np
import yaml
from typing import Dict, Optional
from datetime import datetime


def load_route_prior(config_path: Path) -> Dict:
    """Load route prior configuration."""
    with open(config_path, 'r') as f:
        return yaml.safe_load(f)


def check_sequence_order(events: pd.DataFrame, route_prior: Dict) -> Dict:
    """
    Check if detected landmarks are in correct sequence order.

    Parameters
    ----------
    events : pd.DataFrame
        Final landmark events
    route_prior : Dict
        Route configuration with expected order

    Returns
    -------
    dict
        QC metrics for sequence ordering
    """
    expected_order = [lm['id'] for lm in route_prior['landmarks']]
    detected_order = events.sort_values('timestamp_sec')['landmark_id'].tolist()

    # Check if detected order is a subsequence of expected order
    expected_idx = 0
    correct_order = True
    order_errors = []

    for det_lm in detected_order:
        if det_lm in expected_order[expected_idx:]:
            # Find position in expected order
            new_idx = expected_order[expected_idx:].index(det_lm) + expected_idx
            expected_idx = new_idx + 1
        else:
            correct_order = False
            order_errors.append(f"{det_lm} out of order")

    return {
        'correct_order': correct_order,
        'expected_order': expected_order,
        'detected_order': detected_order,
        'order_errors': order_errors,
        'sequence_accuracy': 1.0 if correct_order else 0.0
    }


def compute_detection_rate(events: pd.DataFrame, route_prior: Dict) -> Dict:
    """
    Compute detection rate: fraction of landmarks found.

    Parameters
    ----------
    events : pd.DataFrame
        Final landmark events
    route_prior : Dict
        Route configuration

    Returns
    -------
    dict
        Detection rate metrics
    """
    expected_landmarks = set([lm['id'] for lm in route_prior['landmarks']])
    detected_landmarks = set(events['landmark_id'].unique())

    detected_count = len(detected_landmarks)
    expected_count = len(expected_landmarks)
    detection_rate = detected_count / expected_count if expected_count > 0 else 0.0

    missing = expected_landmarks - detected_landmarks
    spurious = detected_landmarks - expected_landmarks

    return {
        'detected_count': detected_count,
        'expected_count': expected_count,
        'detection_rate': detection_rate,
        'missing_landmarks': sorted(list(missing)),
        'spurious_landmarks': sorted(list(spurious))
    }


def compute_arc_length_alignment(events: pd.DataFrame, route_prior: Dict) -> Dict:
    """
    Evaluate arc-length alignment quality.

    Parameters
    ----------
    events : pd.DataFrame
        Final landmark events with s_m column
    route_prior : Dict
        Route configuration with expected arc-lengths

    Returns
    -------
    dict
        Arc-length alignment metrics
    """
    lm_lookup = {lm['id']: lm for lm in route_prior['landmarks']}

    errors = []
    for _, event in events.iterrows():
        lm_id = event['landmark_id']
        if lm_id not in lm_lookup:
            continue

        s_expected = lm_lookup[lm_id]['s_m']
        s_detected = event['s_m']
        s_std = lm_lookup[lm_id]['s_std']

        error = abs(s_detected - s_expected)
        normalized_error = error / s_std

        errors.append({
            'landmark_id': lm_id,
            'error_m': error,
            'normalized_error': normalized_error
        })

    if not errors:
        return {'mean_error_m': 0, 'mean_normalized_error': 0, 'errors': []}

    df_errors = pd.DataFrame(errors)

    return {
        'mean_error_m': df_errors['error_m'].mean(),
        'std_error_m': df_errors['error_m'].std(),
        'mean_normalized_error': df_errors['normalized_error'].mean(),
        'max_error_m': df_errors['error_m'].max(),
        'errors': errors
    }


def generate_report(events: pd.DataFrame, route_prior: Dict, output_path: Path):
    """
    Generate QC report markdown file.

    Parameters
    ----------
    events : pd.DataFrame
        Final landmark events
    route_prior : Dict
        Route configuration
    output_path : Path
        Output path for report
    """
    # Compute metrics
    sequence_metrics = check_sequence_order(events, route_prior)
    detection_metrics = compute_detection_rate(events, route_prior)
    alignment_metrics = compute_arc_length_alignment(events, route_prior)

    # Generate report
    report_lines = [
        "# Landmark Sequence QC Report",
        f"\n**Generated:** {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
        f"\n**Route:** {route_prior['route']['name']}",
        "\n---\n",
        "## Summary\n",
        f"- **Total landmarks detected:** {len(events)}",
        f"- **Detection rate:** {detection_metrics['detection_rate']:.2%} "
        f"({detection_metrics['detected_count']}/{detection_metrics['expected_count']})",
        f"- **Sequence accuracy:** {sequence_metrics['sequence_accuracy']:.2%}",
        f"- **Mean arc-length error:** {alignment_metrics.get('mean_error_m', 0):.2f} m",
        "\n---\n",
        "## Detection Rate\n",
    ]

    if detection_metrics['missing_landmarks']:
        report_lines.append(f"**Missing landmarks:** {', '.join(detection_metrics['missing_landmarks'])}")
    else:
        report_lines.append("**All expected landmarks detected!**")

    if detection_metrics['spurious_landmarks']:
        report_lines.append(f"\n**Spurious detections:** {', '.join(detection_metrics['spurious_landmarks'])}")

    report_lines.extend([
        "\n---\n",
        "## Sequence Order\n",
        f"**Expected order:** {' → '.join(sequence_metrics['expected_order'])}",
        f"\n**Detected order:** {' → '.join(sequence_metrics['detected_order'])}",
    ])

    if sequence_metrics['correct_order']:
        report_lines.append("\n✅ **Sequence order is correct!**")
    else:
        report_lines.append(f"\n❌ **Sequence order errors:**")
        for error in sequence_metrics['order_errors']:
            report_lines.append(f"- {error}")

    report_lines.extend([
        "\n---\n",
        "## Arc-Length Alignment\n",
        f"- **Mean error:** {alignment_metrics.get('mean_error_m', 0):.2f} ± "
        f"{alignment_metrics.get('std_error_m', 0):.2f} m",
        f"- **Max error:** {alignment_metrics.get('max_error_m', 0):.2f} m",
        f"- **Mean normalized error:** {alignment_metrics.get('mean_normalized_error', 0):.2f} σ",
        "\n### Per-Landmark Errors\n",
    ])

    if alignment_metrics['errors']:
        report_lines.append("| Landmark | Error (m) | Normalized Error (σ) |")
        report_lines.append("|----------|-----------|----------------------|")
        for err in alignment_metrics['errors']:
            report_lines.append(
                f"| {err['landmark_id']} | {err['error_m']:.2f} | {err['normalized_error']:.2f} |"
            )

    report_lines.extend([
        "\n---\n",
        "## Detected Events\n",
        "| Landmark | Timestamp (s) | Arc-length (m) | Score | Method |",
        "|----------|---------------|----------------|-------|--------|"
    ])

    for _, event in events.iterrows():
        report_lines.append(
            f"| {event['landmark_id']} | {event['timestamp_sec']:.2f} | "
            f"{event['s_m']:.1f} | {event['score']:.3f} | {event['method']} |"
        )

    report_lines.extend([
        "\n---\n",
        "## QA Gates\n",
        f"- **Sequence accuracy ≥ 0.80:** {'✅ PASS' if sequence_metrics['sequence_accuracy'] >= 0.80 else '❌ FAIL'}",
        f"- **Detection rate ≥ 0.70:** {'✅ PASS' if detection_metrics['detection_rate'] >= 0.70 else '❌ FAIL'}",
    ])

    # Write report
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, 'w') as f:
        f.write('\n'.join(report_lines))

    # Print summary to console
    print("\n" + "="*70)
    print("QC Report Summary")
    print("="*70)
    print(f"Detection rate: {detection_metrics['detection_rate']:.2%}")
    print(f"Sequence accuracy: {sequence_metrics['sequence_accuracy']:.2%}")
    print(f"Mean arc-length error: {alignment_metrics.get('mean_error_m', 0):.2f} m")
    print("="*70)


def main():
    """Main execution."""
    print("=" * 70)
    print("QC Report Generation")
    print("=" * 70)

    # Load configuration
    route_prior = load_route_prior(Path('conf/route_prior.yaml'))

    # Load events
    events_path = Path('out/events.csv')
    if not events_path.exists():
        print(f"Error: Events file not found at {events_path}")
        print("Run fuse_and_emit_events.py first")
        return

    events = pd.read_csv(events_path)
    print(f"\nLoaded {len(events)} landmark events")

    # Generate report
    output_path = Path('out/qc/sequence_qc_report.md')
    generate_report(events, route_prior, output_path)

    print(f"\n✓ Saved QC report to: {output_path}")


if __name__ == "__main__":
    main()
