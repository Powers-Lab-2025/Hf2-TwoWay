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

    u_latest = sn.Universe(latest)
    u_prev = sn.Universe(previous)

    mols_latest = u_latest.molecules
    mols_prev = u_prev.molecules

    # Check for failure
    if frame_counter % ORIENT_FAIL_INTERVAL == 0:
        if check_failure(mols_latest):
            return 2, "Failure: S < 0.6", ""

    # Check for hops/spinoffs
    if frame_counter % HOP_CHECK_INTERVAL == 0:
        if frame_counter - last_global_spin_frame < COOLDOWN_GLOBAL:
            return 4, "", ""

        spinoff_targets = check_hop(mols_latest, mols_prev)

        for mol_id, dist, from_col in spinoff_targets:
            last_spin = molecule_last_spin_frame.get(mol_id, -COOLDOWN_MOLECULE)
            if frame_counter - last_spin >= COOLDOWN_MOLECULE:
                frame_num = extract_frame_number(latest.name)
                filename = f"{REF_XYZ_PREFIX}_m{mol_id}_t{frame_num}.dyn"
                logline = f"Spinoff: Mol {mol_id} hopped from column {from_col} at frame {frame_num} (dist={dist:.2f} A)"

                last_global_spin_frame = frame_counter
                molecule_last_spin_frame[mol_id] = frame_counter

                return 1, logline, filename

    return 4, "", ""

def check_failure(mols):
    return sn.compute_orientational_order_parameter(mols) < 0.6

def check_hop(mols_latest, mols_prev):
    coms_prev = np.array([mol.center_of_mass() for mol in mols_prev])
    coms_latest = np.array([mol.center_of_mass() for mol in mols_latest])
    box = mols_latest[0].box

    anchors_prev, labels_prev = sn.guess_column_positions(coms_prev, box, return_noise=True)
    assigned_columns, r = sn.assign_fragments_to_columns(coms_prev, anchors_prev, box)
    anchors_latest, labels_latest = sn.guess_column_positions(coms_latest, box, return_noise=True)

    spin_offs = []
    for i, label in enumerate(labels_latest):
        if label == -1:
            anchor = anchors_prev[assigned_columns[i]]
            box_xy = box[:2]
            dist = np.linalg.norm((coms_latest[i][:2] - anchor + box_xy / 2) % box_xy - box_xy / 2)
            if dist > HOP_DISTANCE_FACTOR * r:
                mol_id = mols_latest[i].id
                spin_offs.append((mol_id, dist, assigned_columns[i]))

    spin_offs.sort(key=lambda x: -x[1])
    return spin_offs[:MAX_SPINOFFS_PER_STEP]

def extract_frame_number(filename):
    parts = filename.split(".")
    if len(parts) > 1 and parts[-1].isdigit():
        return f"{int(parts[-1]):05d}"
    return "00000"
