import numpy as np
from pathlib import Path
from sklearn.cluster import DBSCAN
from sn5 import *
from hf2.config import (
    REF_XYZ_PREFIX,
    ORIENT_FAIL_INTERVAL,
    HOP_CHECK_INTERVAL,
    HOP_DISTANCE_FACTOR,
    MAX_SPINOFFS_PER_STEP,
    COOLDOWN_GLOBAL,
    COOLDOWN_MOLECULE,
)

frame_counter = 0
last_global_spin_frame = -COOLDOWN_GLOBAL
molecule_last_spin_frame = {}  # mol_id (1-indexed) -> last frame spun off

def analysis(sim_dir):
    global frame_counter, last_global_spin_frame, molecule_last_spin_frame

    sim_dir = Path(sim_dir)
    traj_files = sorted([
        f for f in sim_dir.iterdir()
        if f.is_file()
        and f.name.startswith(f"{REF_XYZ_PREFIX}.")
        and len(f.suffixes) == 1
        and f.suffix[1:].isdigit()
    ], key=lambda f: f.stat().st_mtime)

    if len(traj_files) < 2:
        return 4, "", ""

    latest = traj_files[-1]
    previous = traj_files[-2]
    frame_counter += 1

    read_frame(str(latest))
    mol_ids = np.array(MOL)  # same length as COORDS.T
    cc_latest = aa_cc()
    box = BOXXYZ[:3]

    # Orientational order failure check
    if frame_counter % ORIENT_FAIL_INTERVAL == 0:
        S = oriS(verbose=False)
        if S < 0.6:
            return 2, "Failure: S < 0.6", ""

    # Hop detection check
    if frame_counter % HOP_CHECK_INTERVAL == 0:
        if frame_counter - last_global_spin_frame < COOLDOWN_GLOBAL:
            return 4, "", ""

        read_frame(str(previous))
        cc_prev = aa_cc()
        anchors_prev, labels_prev = guess_column_positions(cc_prev, box)
        assigned_columns, r = assign_molecules_to_columns(cc_prev, anchors_prev, box)

        read_frame(str(latest))
        cc_latest = aa_cc()
        anchors_latest, labels_latest = guess_column_positions(cc_latest, box)

        spin_offs = []
        for i, label in enumerate(labels_latest):
            if label == -1:
                anchor = anchors_prev[assigned_columns[i]]
                box_xy = box[:2]
                dist = np.linalg.norm((cc_latest[:, i][:2] - anchor + box_xy / 2) % box_xy - box_xy / 2)
                if dist > HOP_DISTANCE_FACTOR * r:
                    spin_offs.append((i, dist, assigned_columns[i]))

        spin_offs.sort(key=lambda x: -x[1])
        spin_offs = spin_offs[:MAX_SPINOFFS_PER_STEP]

        for mol_idx, dist, from_col in spin_offs:
            mol_id = mol_ids[mol_idx]  # 1-indexed
            last_spin = molecule_last_spin_frame.get(mol_id, -COOLDOWN_MOLECULE)
            if frame_counter - last_spin >= COOLDOWN_MOLECULE:
                frame_num = extract_frame_number(latest.name)
                filename = f"{REF_XYZ_PREFIX}_m{mol_id}_t{frame_num}.dyn"
                logline = f"Spinoff: Mol {mol_id} (idx {mol_idx}) hopped from column {from_col} at frame {frame_num} (dist={dist:.2f} A)"
                last_global_spin_frame = frame_counter
                molecule_last_spin_frame[mol_id] = frame_counter
                return 1, logline, filename

    return 4, "", ""

def extract_frame_number(filename):
    parts = filename.split(".")
    if len(parts) > 1 and parts[-1].isdigit():
        return f"{int(parts[-1]):05d}"
    return "00000"

def guess_column_positions(coms, box, eps=2, min_samples=5):
    xy_coords = coms[:2, :].T % box[:2]
    db = DBSCAN(eps=eps, min_samples=min_samples).fit(xy_coords)
    labels = db.labels_
    anchors = []
    for label in sorted(set(labels)):
        if label == -1:
            continue
        cluster = xy_coords[labels == label]
        anchors.append(np.mean(cluster, axis=0))
    return np.array(anchors), labels

def assign_molecules_to_columns(coms, anchors, box):
    xy_coms = coms[:2, :].T % box[:2]
    assignments = []
    for xy in xy_coms:
        dists = np.linalg.norm((anchors - xy + box[:2]/2) % box[:2] - box[:2]/2, axis=1)
        assignments.append(np.argmin(dists))
    assignments = np.array(assignments)
    # Estimate r as average nearest-neighbor column spacing
    dists = [
        np.min([
            np.linalg.norm((a - b + box[:2]/2) % box[:2] - box[:2]/2)
            for j, b in enumerate(anchors) if i != j
        ])
        for i, a in enumerate(anchors)
    ]
    r = np.mean(dists)
    return assignments, r
