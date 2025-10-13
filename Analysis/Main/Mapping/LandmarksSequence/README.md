# Landmark Sequence Recognition with Route Priors

This pipeline adds an **offline re-ranking layer** that turns noisy vision detections into **ordered landmark events** by leveraging:
1. Fixed walk order from route prior
2. Dead-reckoned arc-length along the route
3. HMM-based sequence alignment

## Pipeline Overview

```
Vision Detections → Arc-Length → Beam Filter → HMM Align → Fuse → Events
                                     ↓              ↓
                                events_beam    events_hmm
```

## Directory Structure

```
LandmarksSequence/
  conf/
    route_prior.yaml          # Expected landmark order and arc-lengths
    sequence_priors.yaml      # Detection gates and fusion weights
  scripts/
    build_arc_length.py       # Compute cumulative distance from XSens
    ordered_beam_filter.py    # Greedy order-based filtering
    hmm_sequence_align.py     # HMM/Viterbi sequence decoding
    text_features.py          # OCR/text analysis for signs, boards, plaques
    fuse_and_emit_events.py   # Merge beam + HMM results
    qc_report.py              # Evaluate sequence accuracy
  out/
    arc_length.csv            # Cumulative distance ŝ(t)
    events_beam.csv           # Beam filter results
    events_hmm.csv            # HMM results
    events.csv                # Final fused events
    qc/
      sequence_qc_report.md   # Quality control report
```

## Usage

### 1. Build Arc-Length Function

Compute cumulative distance from XSens COM trajectory:

```bash
cd Analysis/Main/Mapping/LandmarksSequence
python scripts/build_arc_length.py
```

**Input:** `../data/xsens_com_enu.csv`
**Output:** `out/arc_length.csv`

### 2. Run Ordered Beam Filter

Greedy baseline using order constraints and distance windows:

```bash
python scripts/ordered_beam_filter.py
```

**Input:**
- `conf/route_prior.yaml`
- `conf/sequence_priors.yaml`
- `../Landmarks/out/*.csv` (existing detection results)
- `out/arc_length.csv`

**Output:** `out/events_beam.csv`

### 3. Run HMM Sequence Alignment

HMM/Viterbi decoding for robust sequence alignment:

```bash
python scripts/hmm_sequence_align.py
```

**Input:** Same as beam filter
**Output:** `out/events_hmm.csv`

### 4. Fuse Results

Merge beam and HMM results to emit final events:

```bash
python scripts/fuse_and_emit_events.py
```

**Input:** `out/events_beam.csv`, `out/events_hmm.csv`
**Output:** `out/events.csv`

### 5. Generate QC Report

Evaluate sequence accuracy and detection quality:

```bash
python scripts/qc_report.py
```

**Input:** `out/events.csv`, `conf/route_prior.yaml`
**Output:** `out/qc/sequence_qc_report.md`

## Configuration

### Route Prior (`conf/route_prior.yaml`)

Defines expected landmark order and arc-length positions:

```yaml
route:
  name: "UCLA Walk v1"
  start_desc: "NPI lobby doors"

landmarks:
  - id: LM001
    name: "History display case"
    s_m: 65.0        # Expected arc-length position (meters)
    s_std: 30.0      # Uncertainty (meters)
    category: lookalike
```

### Sequence Priors (`conf/sequence_priors.yaml`)

Detection gates and fusion weights:

```yaml
gates:
  s_window_m: 100.0              # Arc-length search window
  min_persistence_frames: 12     # Minimum detection duration
  min_score: 0.28                # CLIP score threshold

weights:
  w_vision: 0.50    # CLIP similarity weight
  w_sdist: 0.30     # Arc-length distance weight
  w_ocr: 0.10       # OCR feature weight (future)
  w_layout: 0.05    # Layout fingerprint weight (future)
  w_heading: 0.05   # Heading alignment weight (future)

hmm:
  self_trans: 0.92  # P(stay in same landmark)
  fwd_trans: 0.07   # P(move to next landmark)
  skip_trans: 0.01  # P(skip one landmark)
```

## QA Gates

The pipeline passes QA if:
- ✅ Sequence accuracy ≥ 0.80
- ✅ Look-alike F1 ≥ 0.70
- ✅ 90% of accepted events persist ≥ 12 frames

## Integration with Main Pipeline

The final `out/events.csv` can be integrated with:
1. GPS anchoring for world-frame localization
2. SLAM refinement for loop closure
3. Waypoint PDR for drift correction

## Dependencies

- pandas
- numpy
- scipy
- pyyaml

All dependencies should be available in the main `venv` environment.
