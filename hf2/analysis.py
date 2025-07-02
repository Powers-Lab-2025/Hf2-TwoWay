import numpy as np
from pathlib import Path
import sn5 as sn
from hf2.config import (
    REF_XYZ_PREFIX,
    ORIENT_FAIL_INTERVAL,
    HOP_CHECK_INTERVAL,
    HOP_DISTANCE_FACTOR,
    HOP_REFERENCE_FRAMES_BACK,
    MAX_SPINOFFS_PER_STEP,
    COOLDOWN_GLOBAL,
    COOLDOWN_MOLECULE,
)

from sklearn.cluster import DBSCAN

def guess_column_positions(coms, box, n_molecules=600, n_columns=30, return_noise=False, verbose=False):
    """
    Cluster XY COMs to find column centers using DBSCAN.
    
    Args:
        coms: Center of mass coordinates
        box: Simulation box dimensions
        n_molecules: Expected number of molecules (for parameter tuning)
        n_columns: Expected number of columns (for parameter tuning)
        return_noise: Whether to return noise labels
        verbose: Print diagnostic information
    
    Returns: (anchors, labels)
    """
    xy_coords = coms[:2].T  # shape: (N, 2)
    box_xy = box[:2]
    xy_wrapped = xy_coords % box_xy

    # Estimating column spacing from box size and number of columns
    expected_column_spacing = min(box_xy) / np.sqrt(n_columns)
    eps = expected_column_spacing * 0.3
    molecules_per_column = n_molecules // n_columns
    min_samples = max(5, molecules_per_column - 5)

    if verbose:
        print(f"DBSCAN parameters: eps={eps:.2f}, min_samples={min_samples}")
        print(f"Expected column spacing: {expected_column_spacing:.2f}")

    db = DBSCAN(eps=eps, min_samples=min_samples).fit(xy_wrapped)
    labels = db.labels_

    if verbose:
        n_clusters = len(set(labels)) - (1 if -1 in labels else 0)
        n_noise = np.sum(labels == -1)
        print(f"Found {n_clusters} columns, {n_noise} noise molecules")
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
        print(f"[ANALYSIS] Average column spacing r = {r:.2f} A")

    return assigned_columns, r


def analysis(sim_dir, path_state):
    """
    Analysis function that uses path_state for per-path tracking instead of globals.
    
    Args:
        sim_dir: Path to simulation directory
        path_state: SimulationPath object containing state variables
    
    Returns:
        tuple: (code, logline, filename)
    """
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

    if len(traj_files) < HOP_REFERENCE_FRAMES_BACK + 1:
        return 4, "", ""

    latest = traj_files[-1]
    reference = traj_files[-(HOP_REFERENCE_FRAMES_BACK + 1)]  # N frames back
    path_state.analysis_frame_counter += 1

    print(f"[ANALYSIS] Latest XYZ file: {latest}.")

    # Read latest frame
    sn.read_frame(latest)
    com_latest = sn.aa_cm()
    box = sn.DIMS.copy()

    if path_state.analysis_frame_counter % ORIENT_FAIL_INTERVAL == 0:
        if sn.oriS() < 0.6:
            return 2, f"Failure: S < 0.6 at frame {extract_frame_number(latest.name)}", ""

    if path_state.analysis_frame_counter % HOP_CHECK_INTERVAL == 0:
        if path_state.analysis_frame_counter - path_state.last_global_spin_frame < COOLDOWN_GLOBAL:
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

        # Check for molecular hops
        for i, label in enumerate(labels_latest):
            if label == -1:
                original_column_idx = assigned_columns[i]
                original_anchor = anchors_ref[original_column_idx]

                box_xy = box[:2]
                dist = np.linalg.norm((com_latest[:2, i] - original_anchor + box_xy / 2) % box_xy - box_xy / 2)

                if dist > HOP_DISTANCE_FACTOR * r:
                    mol_id = i
                    last_spin = path_state.molecule_last_spin_frame.get(mol_id, -COOLDOWN_MOLECULE)

                    if path_state.analysis_frame_counter - last_spin >= COOLDOWN_MOLECULE:
                        frame_num = extract_frame_number(latest.name)
                        filename = f"{REF_XYZ_PREFIX}_m{mol_id}_t{frame_num}.dyn"
                        
                        target_anchor = None
                        target_column_idx = None
                        if len(anchors_latest) > 0:
                            dists_to_current = np.linalg.norm(
                                (anchors_latest - com_latest[:2, i] + box_xy / 2) % box_xy - box_xy / 2, 
                                axis=1
                            )
                            target_column_idx = np.argmin(dists_to_current)
                            target_anchor = anchors_latest[target_column_idx]
                        
                        original_pos_str = f"({original_anchor[0]:.2f}, {original_anchor[1]:.2f})"
                        if target_anchor is not None:
                            target_pos_str = f"({target_anchor[0]:.2f}, {target_anchor[1]:.2f})"
                            logline = f"Spinoff: Mol {mol_id} hopped from column {original_column_idx} {original_pos_str} to near column {target_column_idx} {target_pos_str} at frame {frame_num} (dist={dist:.2f} Å)"
                        else:
                            logline = f"Spinoff: Mol {mol_id} hopped from column {original_column_idx} {original_pos_str} to unassigned position at frame {frame_num} (dist={dist:.2f} Å)"

                        path_state.last_global_spin_frame = path_state.analysis_frame_counter
                        path_state.molecule_last_spin_frame[mol_id] = path_state.analysis_frame_counter
                        return 1, logline, filename

        # If no valid spinoffs
        return 4, "", ""

    return 4, "", ""

def extract_frame_number(filename):
    parts = filename.split(".")
    if len(parts) > 1 and parts[-1].isdigit():
        return f"{int(parts[-1]):05d}"
    return "00000"