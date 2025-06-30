import numpy as np
from pathlib import Path
import MDAnalysis as mda
import ar5 as ar
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
last_global_spin_frame = -COOLDOWN_GLOBAL  # initially allow spin
fragment_last_spin_frame = {}  # frag_idx -> last frame spun off


def analysis(sim_dir):
    """
    Main analysis. Returns a tuple:
    (action_code, log_message, spin_off_filename)

    Codes:
    1 - Spin off
    2 - Stop as failure
    3 - Stop as success
    4 - Continue
    """
    global frame_counter, last_global_spin_frame, fragment_last_spin_frame

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

    u_latest = mda.Universe(str(latest), format="TXYZ")
    u_prev = mda.Universe(str(previous), format="TXYZ")

    frags_latest = list(u_latest.select_atoms("all").fragments)
    frags_prev = list(u_prev.select_atoms("all").fragments)

    # Check for failure
    if frame_counter % ORIENT_FAIL_INTERVAL == 0:
        if check_failure(frags_latest):
            return 2, "Failure: S < 0.6", ""

    # Check for hopsspinoffs
    if frame_counter % HOP_CHECK_INTERVAL == 0:
        # Enforce global cooldown
        if frame_counter - last_global_spin_frame < COOLDOWN_GLOBAL:
            return 4, "", ""

        spinoff_targets = check_hop(frags_latest, frags_prev)

        for frag_idx, dist, from_col in spinoff_targets:
            last_spin = fragment_last_spin_frame.get(frag_idx, -COOLDOWN_MOLECULE)
            if frame_counter - last_spin >= COOLDOWN_MOLECULE:
                mol_id = frags_latest[frag_idx][0].id
                frame_num = extract_frame_number(latest.name)
                filename = f"{REF_XYZ_PREFIX}_m{mol_id}_t{frame_num}.dyn"
                logline = f"Spinoff: Fragment {frag_idx} (Mol {mol_id}) hopped from column {from_col} at frame {frame_num} (dist={dist:.2f} A)"

                # Update cooldowns
                last_global_spin_frame = frame_counter
                fragment_last_spin_frame[frag_idx] = frame_counter

                return 1, logline, filename

    return 4, "", ""


def check_failure(frags):
    S = ar.compute_orientational_order_parameter(frags, verbose=False)
    return S < 0.6


def check_hop(frags_latest, frags_prev):
    coms_prev = np.array([frag.center_of_mass() for frag in frags_prev])
    coms_latest = np.array([frag.center_of_mass() for frag in frags_latest])
    box = frags_latest[0].dimensions[:3]

    anchors_prev, labels_prev = ar.guess_column_positions(coms_prev, box, return_noise=True)
    assigned_columns, r = ar.assign_fragments_to_columns(coms_prev, anchors_prev, box)

    anchors_latest, labels_latest = ar.guess_column_positions(coms_latest, box, return_noise=True)

    spin_offs = []
    for i, label in enumerate(labels_latest):
        if label == -1:
            anchor = anchors_prev[assigned_columns[i]]
            box_xy = box[:2]
            dist = np.linalg.norm((coms_latest[i][:2] - anchor + box_xy / 2) % box_xy - box_xy / 2)
            if dist > HOP_DISTANCE_FACTOR * r:
                spin_offs.append((i, dist, assigned_columns[i]))

    spin_offs.sort(key=lambda x: -x[1])
    return spin_offs[:MAX_SPINOFFS_PER_STEP]


def extract_frame_number(filename):
    parts = filename.split(".")
    if len(parts) > 1 and parts[-1].isdigit():
        return f"{int(parts[-1]):05d}"
    return "00000"
