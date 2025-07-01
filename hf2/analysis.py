import numpy as np
from pathlib import Path
import sn5 as sn
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
molecule_last_spin_frame = {}  # mol_id -> last spin frame


from sklearn.cluster import DBSCAN

def guess_column_positions(coms, box, eps=2.0, min_samples=5, return_noise=False, verbose=False):
    """
    Cluster XY COMs to find column centers.
    Returns: (anchors, labels)
    """
    xy_coords = coms[:2].T  # shape: (N, 2)
    box_xy = box[:2]
    xy_wrapped = xy_coords % box_xy

    db = DBSCAN(eps=eps, min_samples=min_samples).fit(xy_wrapped)
    labels = db.labels_

    if verbose:
        n_noise = np.sum(labels == -1)
        if n_noise > 0:
            print(f"{n_noise} out of {len(labels)} COMs were not assigned to any column.")
        else:
            print("All COMs assigned to columns.")

    anchors = []
    for label in sorted(set(labels)):
        if label == -1:
            continue
        cluster_points = xy_wrapped[labels == label]
        anchor = np.mean(cluster_points, axis=0)
        anchors.append(anchor)

    anchors = np.array(anchors)

    if return_noise:
        return anchors, labels
    else:
        return anchors


def assign_fragments_to_columns(coms, anchors, box, verbose=False):
    """
    Assign each COM to its nearest anchor using periodic XY distance.
    Returns (assigned_columns, r) where r is average inter-anchor spacing.
    """
    xy_coords = coms[:2].T  # shape: (N, 2)
    box_xy = box[:2]
    xy_wrapped = xy_coords % box_xy

    assigned_columns = []
    for xy in xy_wrapped:
        dists = np.linalg.norm((anchors - xy + box_xy / 2) % box_xy - box_xy / 2, axis=1)
        assigned_columns.append(np.argmin(dists))
    assigned_columns = np.array(assigned_columns)

    # Estimate average inter-column spacing r
    nearest_dists = []
    for i, a in enumerate(anchors):
        dists = np.linalg.norm((anchors - a + box_xy / 2) % box_xy - box_xy / 2, axis=1)
        dists[i] = np.inf
        nearest_dists.append(np.min(dists))
    r = np.mean(nearest_dists)

    if verbose:
        print(f"Average column spacing r = {r:.2f} A")

    return assigned_columns, r


def analysis(sim_dir):
    global frame_counter, last_global_spin_frame, molecule_last_spin_frame

    sim_dir = Path(sim_dir)
    traj_files = sorted([
        f for f in sim_dir.iterdir()
        if f.is_file()
        and f.name.startswith(f"{REF_XYZ_PREFIX}.")
        and len(f.suffixes) == 1
        and f.suffix[1:].isdigit()
        and not f.name.endswith(".dyn")
        and not f.name.endswith(".key")
    ], key=lambda f: f.stat().st_mtime)

    if len(traj_files) < 2:
        return 4, "", ""

    latest = traj_files[-1]
    previous = traj_files[-2]
    frame_counter += 1

    # Read latest frame
    sn.read_frame(latest)
    com_latest = sn.aa_cm()
    box = sn.DIMS.copy()

    if frame_counter % ORIENT_FAIL_INTERVAL == 0:
        if sn.oriS() < 0.6:
            return 2, f"Failure: S < 0.6 at frame {extract_frame_number(latest.name)}", ""

    if frame_counter % HOP_CHECK_INTERVAL == 0:
        if frame_counter - last_global_spin_frame < COOLDOWN_GLOBAL:
            return 4, "", ""

        sn.read_frame(previous)
        com_prev = sn.aa_cm()

        # cluster columns on previous frame
        anchors_prev, labels_prev = guess_column_positions(com_prev, box, return_noise=True)
        assigned_columns, r = assign_fragments_to_columns(com_prev, anchors_prev, box)


        # cluster latest frame
        sn.read_frame(latest)
        com_latest = sn.aa_cm()
        anchors_latest, labels_latest = guess_column_positions(com_latest, box, return_noise=True)

        spin_offs = []
        for i, label in enumerate(labels_latest):
            if label == -1:
                anchor = anchors_prev[assigned_columns[i]]
                box_xy = box[:2]
                dist = np.linalg.norm((com_latest[:2, i] - anchor + box_xy / 2) % box_xy - box_xy / 2)

                if dist > HOP_DISTANCE_FACTOR * r:
                    mol_id = i  # index is ID
                    last_spin = molecule_last_spin_frame.get(mol_id, -COOLDOWN_MOLECULE)

                    if frame_counter - last_spin >= COOLDOWN_MOLECULE:
                        frame_num = extract_frame_number(latest.name)
                        filename = f"{REF_XYZ_PREFIX}_m{mol_id}_t{frame_num}.dyn"
                        logline = f"Spinoff: Mol {mol_id} hopped from column {assigned_columns[i]} at frame {frame_num} (dist={dist:.2f} Å)"

                        last_global_spin_frame = frame_counter
                        molecule_last_spin_frame[mol_id] = frame_counter
                        return 1, logline, filename

        # If no valid spinoffs
        return 4, "", ""

    return 4, "", ""

def extract_frame_number(filename):
    parts = filename.split(".")
    if len(parts) > 1 and parts[-1].isdigit():
        return f"{int(parts[-1]):05d}"
    return "00000"
