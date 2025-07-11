import numpy as np
from pathlib import Path
import sn5 as sn
from hf2tw.config import (
    REF_XYZ_PREFIX,
    ORIENT_FAIL_INTERVAL,
)


def pbc_distance(pos1, pos2, box_dims):
    delta = pos1 - pos2
    delta = delta - box_dims * np.round(delta / box_dims)
    return np.linalg.norm(delta)


def pbc_displacement(pos1, pos2, box_dims):
    delta = pos1 - pos2
    delta = delta - box_dims * np.round(delta / box_dims)
    return delta


from scipy.stats import halfnorm


def simple_column_clustering(coms, box, n_columns=30, threshold=0.1, max_iter=100, verbose=False):
    xy_coords = coms[:2].T
    box_xy = box[:2]
    n_molecules = xy_coords.shape[0]
    
    grid_size_x = int(np.ceil(np.sqrt(n_columns * box_xy[0] / box_xy[1])))
    grid_size_y = int(np.ceil(n_columns / grid_size_x))
    
    x_spacing = box_xy[0] / grid_size_x
    y_spacing = box_xy[1] / grid_size_y
    
    anchors = []
    for i in range(grid_size_x):
        for j in range(grid_size_y):
            if len(anchors) < n_columns:
                x = (i + 0.5) * x_spacing - box_xy[0]/2
                y = (j + 0.5) * y_spacing - box_xy[1]/2
                anchors.append([x, y])
    
    anchors = np.array(anchors[:n_columns])
    
    if verbose:
        print(f"[CLUSTERING] Starting with {n_columns} columns")
    
    for iteration in range(max_iter):
        assignments = np.zeros(n_molecules, dtype=int)
        for i, mol_pos in enumerate(xy_coords):
            min_dist = float('inf')
            for j, anchor in enumerate(anchors):
                dist = pbc_distance(mol_pos, anchor, box_xy)
                if dist < min_dist:
                    min_dist = dist
                    assignments[i] = j
        
        new_anchors = np.zeros_like(anchors)
        column_movements = np.zeros(len(anchors))
        
        for col_idx in range(len(anchors)):
            mol_indices = np.where(assignments == col_idx)[0]
            
            if len(mol_indices) == 0:
                new_anchors[col_idx] = anchors[col_idx]
                column_movements[col_idx] = 0
            else:
                ref_pos = xy_coords[mol_indices[0]]
                displacements = []
                for mol_idx in mol_indices:
                    disp = pbc_displacement(xy_coords[mol_idx], ref_pos, box_xy)
                    displacements.append(disp)
                
                displacements = np.array(displacements)
                mean_disp = np.mean(displacements, axis=0)
                new_pos = ref_pos + mean_disp
                
                new_anchors[col_idx] = new_pos
                column_movements[col_idx] = pbc_distance(new_anchors[col_idx], anchors[col_idx], box_xy)
        
        max_movement = np.max(column_movements)
        
        if verbose and iteration % 10 == 0:
            print(f"[CLUSTERING] Iteration {iteration}: max movement = {max_movement:.4f}")
        
        anchors = new_anchors
        
        if max_movement < threshold:
            if verbose:
                print(f"[CLUSTERING] Converged after {iteration + 1} iterations")
            break
    
    labels = assignments.copy()
    
    for col_idx in range(len(anchors)):
        mol_indices = np.where(assignments == col_idx)[0]
        if len(mol_indices) <= 1:
            continue
            
        anchor = anchors[col_idx]
        distances = []
        for mol_idx in mol_indices:
            dist = pbc_distance(xy_coords[mol_idx], anchor, box_xy)
            distances.append(dist)
        
        distances = np.array(distances)
        mean_dist = np.mean(distances)
        std_dist = np.std(distances)
        
        if std_dist > 0:
            for i, mol_idx in enumerate(mol_indices):
                if distances[i] > mean_dist + 3.5 * std_dist:
                    labels[mol_idx] = -1
    
    if verbose:
        n_outliers = np.sum(labels == -1)
        print(f"[CLUSTERING] Found {n_outliers} outliers out of {n_molecules} molecules")
    
    return anchors, labels


def assign_fragments_to_columns(coms, anchors, box, verbose=False):
    xy_coords = coms[:2].T
    box_xy = box[:2]
    
    assigned_columns = []
    for xy in xy_coords:
        min_dist = float('inf')
        best_col = -1
        for j, anchor in enumerate(anchors):
            dist = pbc_distance(xy, anchor, box_xy)
            if dist < min_dist:
                min_dist = dist
                best_col = j
        assigned_columns.append(best_col)
    
    assigned_columns = np.array(assigned_columns)
    
    nearest_dists = []
    for i, anchor_i in enumerate(anchors):
        min_dist = float('inf')
        for j, anchor_j in enumerate(anchors):
            if i != j:
                dist = pbc_distance(anchor_i, anchor_j, box_xy)
                if dist < min_dist:
                    min_dist = dist
        nearest_dists.append(min_dist)
    
    r = np.mean(nearest_dists)
    
    if verbose:
        print(f"[ANALYSIS] Average column spacing r = {r:.2f} A")
    
    return assigned_columns, r


def analysis(sim_dir, path_state):
    sim_dir = Path(sim_dir)
    prefix = path_state.prefix
    
    traj_files = sorted([
        f for f in sim_dir.iterdir()
        if f.is_file()
        and f.name.startswith(f"{prefix}.")
        and len(f.suffixes) == 1
        and f.suffix[1:].isdigit()
        and not f.name.endswith(".dyn")
        and not f.name.endswith(".key")
    ], key=lambda f: f.stat().st_mtime)

    if len(traj_files) < 1:
        return 4, "", ""

    latest = traj_files[-1]
    path_state.analysis_frame_counter += 1

    print(f"[ANALYSIS] Latest XYZ file: {latest}.")

    sn.read_frame(latest)
    com_latest = sn.aa_cm()
    box = sn.DIMS.copy()

    if path_state.analysis_frame_counter % ORIENT_FAIL_INTERVAL == 0:
        if sn.oriS() < 0.6:
            return 2, f"Failure: S < 0.6 at frame {extract_frame_number(latest.name)}", ""

    anchors, labels = simple_column_clustering(com_latest, box, n_columns=30, verbose=False)
    assigned_columns, r = assign_fragments_to_columns(com_latest, anchors, box)
    
    proximity_threshold = 0.1 * r
    
    for mol_id in range(com_latest.shape[1]):
        box_xy = box[:2]
        
        for col_idx, anchor in enumerate(anchors):
            dist = pbc_distance(com_latest[:2, mol_id], anchor, box_xy)
            
            if dist < proximity_threshold:
                if mol_id not in path_state.molecule_displacement_frames:
                    path_state.molecule_displacement_frames[mol_id] = 0
                path_state.molecule_displacement_frames[mol_id] += 1
                
                if path_state.molecule_displacement_frames[mol_id] >= 10:
                    logline = f"Molecule {mol_id+1} within {dist:.2f} A of column {col_idx} anchor for 10 frames. Killing simulation."
                    return 2, logline, ""
                break
        else:
            if mol_id in path_state.molecule_displacement_frames:
                path_state.molecule_displacement_frames[mol_id] = 0

    return 4, "", ""


def extract_frame_number(filename):
    parts = filename.split(".")
    if len(parts) > 1 and parts[-1].isdigit():
        return f"{int(parts[-1]):05d}"
    return "00000"