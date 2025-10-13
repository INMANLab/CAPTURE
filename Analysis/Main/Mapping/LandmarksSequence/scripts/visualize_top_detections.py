#!/usr/bin/env python3
"""
Visualize top detections for each landmark across walks.

For each landmark, shows:
- Original exemplar image(s)
- Top 3 detections from each walk (with scores)

Usage:
    python visualize_top_detections.py --walks Walk1 Walk2 Walk3 Walk4
    python visualize_top_detections.py --all
    python visualize_top_detections.py --landmark LM050  # Single landmark
"""
import sys
from pathlib import Path
import pandas as pd
import numpy as np
import cv2
import yaml
import argparse
from typing import List, Dict, Tuple
import matplotlib.pyplot as plt
import matplotlib.patches as patches

# Add parent to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent))


def load_route_prior(config_path: Path) -> Dict:
    """Load route prior configuration."""
    with open(config_path, 'r') as f:
        return yaml.safe_load(f)


def load_detections_for_walk(detections_dir: Path, walk_name: str) -> pd.DataFrame:
    """Load detection CSV for a specific walk."""
    pattern = f"{walk_name.lower()}_*.csv"
    csv_files = list(detections_dir.glob(pattern))

    if not csv_files:
        return pd.DataFrame()

    return pd.read_csv(csv_files[0])


def find_video_for_walk(walk_name: str, base_dir: Path = None) -> Path:
    """
    Find video file for a walk by searching standard locations.

    Parameters
    ----------
    walk_name : str
        Walk name (e.g., 'Walk1')
    base_dir : Path, optional
        Base directory to search (default: ~/Downloads/RW/RW1)

    Returns
    -------
    Path or None
        Path to video file, or None if not found
    """
    if base_dir is None:
        base_dir = Path.home() / "Downloads/RW/RW1"

    walk_dir = base_dir / walk_name

    if not walk_dir.exists():
        return None

    # Search for Pupil video directories
    pupil_dirs = list(walk_dir.glob("Pupil*"))

    for pupil_dir in pupil_dirs:
        # Find .mp4 files
        video_files = list(pupil_dir.glob("*.mp4"))
        if video_files:
            return video_files[0]  # Return first video found

    return None


def load_video(walk_name: str, base_dir: Path = None) -> Tuple[cv2.VideoCapture, float]:
    """Load video for a walk."""
    video_path = find_video_for_walk(walk_name, base_dir)

    if video_path is None or not video_path.exists():
        return None, 0.0

    cap = cv2.VideoCapture(str(video_path))
    fps = cap.get(cv2.CAP_PROP_FPS)

    return cap, fps


def extract_frame_crop(cap: cv2.VideoCapture, frame_num: int, bbox_str: str) -> np.ndarray:
    """Extract cropped region from video frame."""
    # Parse bbox
    coords = [int(x) for x in bbox_str.split(',')]
    x1, y1, x2, y2 = coords

    # Seek to frame
    cap.set(cv2.CAP_PROP_POS_FRAMES, frame_num)
    ret, frame = cap.read()

    if not ret:
        return None

    # Extract crop
    crop = frame[y1:y2, x1:x2]

    # Convert BGR to RGB
    crop_rgb = cv2.cvtColor(crop, cv2.COLOR_BGR2RGB)

    return crop_rgb


def load_exemplar_images(exemplar_dir: Path, lm_id: str) -> List[np.ndarray]:
    """Load exemplar images for a landmark."""
    lm_dir = exemplar_dir / lm_id

    if not lm_dir.exists():
        return []

    images = []
    for img_path in sorted(lm_dir.glob("*.jpg")) + sorted(lm_dir.glob("*.JPEG")) + sorted(lm_dir.glob("*.png")):
        img = cv2.imread(str(img_path))
        if img is not None:
            img_rgb = cv2.cvtColor(img, cv2.COLOR_BGR2RGB)
            images.append(img_rgb)

    return images


def visualize_landmark_detections(lm_id: str, lm_name: str, walks: List[str],
                                  detections_dir: Path, exemplar_dir: Path,
                                  output_dir: Path, top_n: int = 3, base_dir: Path = None):
    """
    Create visualization for a single landmark across all walks.

    Layout:
    - Top row: Exemplar images
    - Subsequent rows: Top N detections from each walk
    """
    print(f"\nVisualizing {lm_id}: {lm_name}")

    # Load exemplars
    exemplars = load_exemplar_images(exemplar_dir, lm_id)

    if not exemplars:
        print(f"  No exemplars found, skipping")
        return

    # Collect top detections from each walk
    walk_detections = {}

    for walk_name in walks:
        detections = load_detections_for_walk(detections_dir, walk_name)

        if len(detections) == 0:
            continue

        # Filter for this landmark
        lm_dets = detections[detections['landmark_id'] == lm_id].copy()

        if len(lm_dets) == 0:
            continue

        # Get top N by score
        top_dets = lm_dets.nlargest(top_n, 'score')

        # Load video
        cap, fps = load_video(walk_name, base_dir)

        if cap is None:
            print(f"  Warning: Could not load video for {walk_name}")
            continue

        # Extract crops
        crops = []
        for _, det in top_dets.iterrows():
            crop = extract_frame_crop(cap, det['video_frame'], det['bbox_xyxy'])
            if crop is not None:
                crops.append({
                    'image': crop,
                    'score': det['score'],
                    'frame': det['video_frame'],
                    'timestamp': det['timestamp_sec']
                })

        cap.release()

        if crops:
            walk_detections[walk_name] = crops

    if not walk_detections:
        print(f"  No detections found across walks")
        return

    # Create figure
    n_walks = len(walk_detections)
    n_cols = max(len(exemplars), top_n)
    n_rows = 1 + n_walks  # Exemplars + walks

    fig = plt.figure(figsize=(4 * n_cols, 4 * n_rows))

    # Plot exemplars
    for i, exemplar in enumerate(exemplars[:n_cols]):
        ax = plt.subplot(n_rows, n_cols, i + 1)
        ax.imshow(exemplar)
        ax.set_title(f"Exemplar {i+1}", fontsize=10, fontweight='bold')
        ax.axis('off')

    # Plot detections for each walk
    for row_idx, (walk_name, crops) in enumerate(sorted(walk_detections.items())):
        for col_idx, crop_data in enumerate(crops[:n_cols]):
            subplot_idx = (row_idx + 1) * n_cols + col_idx + 1  # +1 because row 0 is exemplars
            ax = plt.subplot(n_rows, n_cols, subplot_idx)
            ax.imshow(crop_data['image'])

            title = (f"{walk_name}\n"
                    f"Score: {crop_data['score']:.3f}\n"
                    f"Frame: {crop_data['frame']}\n"
                    f"Time: {crop_data['timestamp']:.1f}s")

            ax.set_title(title, fontsize=8)
            ax.axis('off')

            # Add border colored by score (green=high, red=low)
            score_color = plt.cm.RdYlGn(crop_data['score'])
            for spine in ax.spines.values():
                spine.set_edgecolor(score_color)
                spine.set_linewidth(3)

    # Overall title
    fig.suptitle(f"{lm_id}: {lm_name}", fontsize=16, fontweight='bold', y=0.98)

    plt.tight_layout()

    # Save
    output_path = output_dir / f"{lm_id}_detections.png"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()

    print(f"  ✓ Saved to {output_path}")


def main():
    """Main execution."""
    parser = argparse.ArgumentParser(description='Visualize top detections for landmarks')
    parser.add_argument('--walks', nargs='+', help='Walk names (e.g., Walk1 Walk2)')
    parser.add_argument('--all', action='store_true', help='Process all walks')
    parser.add_argument('--landmark', type=str, help='Process single landmark (e.g., LM050)')
    parser.add_argument('--top-n', type=int, default=3, help='Number of top detections per walk (default: 3)')
    parser.add_argument('--base-dir', type=str, help='Base directory for video files (default: ~/Downloads/RW/RW1)')
    args = parser.parse_args()

    # Parse base directory
    base_dir = Path(args.base_dir).expanduser() if args.base_dir else None

    print("=" * 70)
    print("Top Detections Visualizer")
    print("=" * 70)

    # Load route prior
    route_prior = load_route_prior(Path('conf/route_prior.yaml'))

    # Determine walks
    detections_dir = Path('../Landmarks/out/detections')

    if args.all:
        detection_files = list(detections_dir.glob('walk*_*.csv'))
        walks = sorted(set([f.name.split('_')[0].capitalize() for f in detection_files]))
    elif args.walks:
        walks = args.walks
    else:
        walks = ['Walk1', 'Walk2', 'Walk3', 'Walk4']

    print(f"Processing walks: {walks}")

    # Determine landmarks
    if args.landmark:
        landmarks = [lm for lm in route_prior['landmarks'] if lm['id'] == args.landmark]
        if not landmarks:
            print(f"Error: Landmark {args.landmark} not found in route_prior.yaml")
            sys.exit(1)
    else:
        landmarks = route_prior['landmarks']

    # Process each landmark
    exemplar_dir = Path('../Landmarks/data/exemplars')
    output_dir = Path('out/visualizations')

    for lm_info in landmarks:
        visualize_landmark_detections(
            lm_info['id'],
            lm_info['name'],
            walks,
            detections_dir,
            exemplar_dir,
            output_dir,
            top_n=args.top_n,
            base_dir=base_dir
        )

    print(f"\n{'='*70}")
    print(f"✓ Visualizations saved to: {output_dir}")
    print(f"{'='*70}")


if __name__ == "__main__":
    main()
