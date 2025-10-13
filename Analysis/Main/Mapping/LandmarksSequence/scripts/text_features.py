#!/usr/bin/env python3
"""
Text-based feature extraction for landmark disambiguation.

Extracts and compares text features from detections to improve landmark recognition,
especially for text-bearing landmarks like:
- Signs (street signs, building signs, directional signs)
- Plaques and information boards
- Door labels and room numbers
- Display cases with text
- Any landmark with visible text content

Features:
1. OCR rare-token matching - Uncommon words are more discriminative
2. Text layout fingerprints - Spatial arrangement of text regions
3. Font/style consistency - Visual text appearance

Input:
    - Video frames with bounding boxes
    - Reference text from exemplars
    - Detection CSV with bbox_xyxy

Output:
    - Text similarity scores per detection
"""
import sys
from pathlib import Path
import cv2
import pandas as pd
import numpy as np
from typing import Dict, List, Tuple
import re
from collections import Counter

# OCR library - using easyocr for simplicity
try:
    import easyocr
    HAS_OCR = True
except ImportError:
    HAS_OCR = False
    print("Warning: easyocr not installed. Text features will be disabled.")
    print("Install with: pip install easyocr")


class TextFeatureExtractor:
    """Extract text-based features for landmark matching."""

    def __init__(self):
        """Initialize OCR reader if available."""
        self.reader = None
        if HAS_OCR:
            # Initialize EasyOCR reader (English only for now)
            self.reader = easyocr.Reader(['en'], gpu=False)

    def extract_text(self, image: np.ndarray) -> List[Tuple[str, float, List]]:
        """
        Extract text from image using OCR.

        Parameters
        ----------
        image : np.ndarray
            Input image (BGR format)

        Returns
        -------
        List[Tuple[str, float, List]]
            List of (text, confidence, bbox) tuples
        """
        if not self.reader:
            return []

        # EasyOCR expects RGB
        image_rgb = cv2.cvtColor(image, cv2.COLOR_BGR2RGB)

        # Run OCR
        results = self.reader.readtext(image_rgb)

        return results

    def compute_rare_tokens(self, text_list: List[str]) -> Dict[str, float]:
        """
        Compute rare token scores based on word frequency.

        Rare/unique words are more discriminative for landmark matching.

        Parameters
        ----------
        text_list : List[str]
            List of text strings

        Returns
        -------
        Dict[str, float]
            Token to rarity score (inverse frequency)
        """
        # Tokenize and normalize
        all_tokens = []
        for text in text_list:
            # Remove punctuation and convert to lowercase
            tokens = re.findall(r'\b\w+\b', text.lower())
            all_tokens.extend(tokens)

        # Count token frequencies
        token_counts = Counter(all_tokens)
        total_tokens = len(all_tokens)

        # Compute rarity scores (inverse frequency)
        rarity_scores = {}
        for token, count in token_counts.items():
            freq = count / total_tokens
            rarity_scores[token] = 1.0 / (freq + 1e-6)

        return rarity_scores

    def text_similarity(self, text1: str, text2: str, rarity_scores: Dict[str, float] = None) -> float:
        """
        Compute text similarity with optional rare-token weighting.

        Parameters
        ----------
        text1, text2 : str
            Texts to compare
        rarity_scores : Dict[str, float], optional
            Token rarity scores for weighting

        Returns
        -------
        float
            Similarity score [0, 1]
        """
        # Tokenize
        tokens1 = set(re.findall(r'\b\w+\b', text1.lower()))
        tokens2 = set(re.findall(r'\b\w+\b', text2.lower()))

        if not tokens1 or not tokens2:
            return 0.0

        # Compute weighted Jaccard similarity
        intersection = tokens1 & tokens2
        union = tokens1 | tokens2

        if rarity_scores:
            # Weight by rarity
            inter_weight = sum(rarity_scores.get(t, 1.0) for t in intersection)
            union_weight = sum(rarity_scores.get(t, 1.0) for t in union)
            similarity = inter_weight / (union_weight + 1e-6)
        else:
            # Standard Jaccard
            similarity = len(intersection) / (len(union) + 1e-6)

        return similarity

    def layout_fingerprint(self, text_results: List[Tuple[str, float, List]]) -> np.ndarray:
        """
        Compute text layout fingerprint from OCR results.

        Captures spatial arrangement of text regions.

        Parameters
        ----------
        text_results : List[Tuple[str, float, List]]
            OCR results (text, confidence, bbox)

        Returns
        -------
        np.ndarray
            Layout fingerprint vector
        """
        if not text_results:
            return np.zeros(16)  # Empty fingerprint

        # Extract bounding box centers and sizes
        centers = []
        areas = []

        for text, conf, bbox in text_results:
            # bbox is list of 4 points [(x1,y1), (x2,y2), (x3,y3), (x4,y4)]
            bbox_array = np.array(bbox)
            center = bbox_array.mean(axis=0)
            area = cv2.contourArea(bbox_array)

            centers.append(center)
            areas.append(area)

        centers = np.array(centers)
        areas = np.array(areas)

        # Normalize centers to [0, 1]
        if len(centers) > 0:
            centers_norm = (centers - centers.min(axis=0)) / (centers.ptp(axis=0) + 1e-6)
        else:
            centers_norm = centers

        # Compute layout statistics
        fingerprint = np.concatenate([
            centers_norm.mean(axis=0) if len(centers_norm) > 0 else np.zeros(2),  # Mean position
            centers_norm.std(axis=0) if len(centers_norm) > 0 else np.zeros(2),   # Spread
            [len(text_results)],  # Number of text regions
            [np.mean(areas)] if len(areas) > 0 else [0],  # Mean area
            [np.std(areas)] if len(areas) > 0 else [0],   # Area variance
            # Spatial distribution features
            centers_norm[:5].flatten() if len(centers_norm) >= 5 else np.zeros(10)  # Top 5 positions
        ])

        # Pad or trim to fixed length
        fingerprint = fingerprint[:16]
        if len(fingerprint) < 16:
            fingerprint = np.pad(fingerprint, (0, 16 - len(fingerprint)))

        return fingerprint

    def layout_similarity(self, fp1: np.ndarray, fp2: np.ndarray) -> float:
        """
        Compute similarity between layout fingerprints.

        Parameters
        ----------
        fp1, fp2 : np.ndarray
            Layout fingerprints

        Returns
        -------
        float
            Similarity score [0, 1]
        """
        # Cosine similarity
        norm1 = np.linalg.norm(fp1)
        norm2 = np.linalg.norm(fp2)

        if norm1 == 0 or norm2 == 0:
            return 0.0

        similarity = np.dot(fp1, fp2) / (norm1 * norm2)
        return max(0.0, similarity)  # Clamp to [0, 1]


def compute_text_scores(detections: pd.DataFrame, video_path: Path,
                       exemplar_text: Dict[str, str] = None) -> pd.DataFrame:
    """
    Compute text similarity scores for all detections.

    Parameters
    ----------
    detections : pd.DataFrame
        Detection CSV with bbox_xyxy column
    video_path : Path
        Path to video file
    exemplar_text : Dict[str, str], optional
        Reference text for each landmark

    Returns
    -------
    pd.DataFrame
        Detections with added text_score column
    """
    if not HAS_OCR:
        print("OCR not available - skipping text features")
        detections['text_score'] = 0.0
        return detections

    extractor = TextFeatureExtractor()
    cap = cv2.VideoCapture(str(video_path))

    text_scores = []

    for idx, row in detections.iterrows():
        frame_num = int(row['video_frame'])
        bbox_str = row['bbox_xyxy']
        x1, y1, x2, y2 = map(int, bbox_str.split(','))

        # Read frame and crop
        cap.set(cv2.CAP_PROP_POS_FRAMES, frame_num)
        ret, frame = cap.read()

        if not ret:
            text_scores.append(0.0)
            continue

        crop = frame[y1:y2, x1:x2]

        # Extract text
        text_results = extractor.extract_text(crop)

        # Compute similarity to exemplar text if available
        if exemplar_text and row['landmark_id'] in exemplar_text:
            detected_text = ' '.join([t[0] for t in text_results])
            ref_text = exemplar_text[row['landmark_id']]
            score = extractor.text_similarity(detected_text, ref_text)
        else:
            # No reference - just use text presence as score
            score = 1.0 if len(text_results) > 0 else 0.0

        text_scores.append(score)

    cap.release()

    detections['text_score'] = text_scores
    return detections


if __name__ == "__main__":
    # Example usage
    print("Text Feature Extractor")
    print("=" * 70)

    if not HAS_OCR:
        print("\nERROR: easyocr not installed")
        print("Install with: pip install easyocr")
        sys.exit(1)

    # Example: Extract text from a test image
    extractor = TextFeatureExtractor()

    # Test text similarity
    text1 = "School of Medicine Entrance"
    text2 = "Medical School Entry"
    text3 = "Random Building"

    print(f"\nText similarity examples:")
    print(f"  '{text1}' vs '{text2}': {extractor.text_similarity(text1, text2):.3f}")
    print(f"  '{text1}' vs '{text3}': {extractor.text_similarity(text1, text3):.3f}")
