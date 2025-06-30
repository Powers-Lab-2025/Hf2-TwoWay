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
)

frame_counter = 0

def analysis(sim_dir):
    """
    Main analysis entry point. Called by a conductor object to determine what action to take.
    It performs several modular checks and returns an action code:
    1 - Spin off
    2 - Stop as failure
    3 - Stop as success
    4 - Continue
    """
    global frame_counter
    sim_dir = Path(sim_dir)
    traj_files = sorted(sim_dir.glob(f"{REF_XYZ_PREFIX}.*"), key=lambda f: f.stat().st_mtime)

    if len(traj_files) < 2:
        return 4  # Not enough frames yet

    latest = traj_files[-1]
    previous = traj_files[-2]
    frame_counter += 1

    u_latest = mda.Universe(str(latest), format="TXYZ")
    u_prev = mda.Universe(str(previous), format="TXYZ")

    frags_latest = list(u_latest.select_atoms("all").fragments)
    frags_prev = list(u_prev.select_atoms("all").fragments)

    # Check for failure based on orientational order
    if frame_counter % ORIENT_FAIL_INTERVAL == 0:
        if check_failure(frags_latest):
            log_action("Failure: S < 0.6")
            return 2

    # Check for success/failure or hop opportunities
    if frame_counter % HOP_CHECK_INTERVAL == 0:
        focus_result = check_focus(sim_dir, frags_latest)
        if focus_result is not None:
            return focus_result

        spinoff_targets = check_hop(frags_latest, frags_prev, sim_dir)
        if spinoff_targets:
            for frag_idx, dist, col in spinoff_targets:
                write_focus_file(sim_dir, frag_idx, col, dist, frame_number=frame_counter)
                log_action(f"Spinoff: Fragment {frag_idx} is {dist:.2f} A from column {col}")
            return 1

    return 4


def check_failure(frags):
    """Returns True if orientational order parameter S < 0.6, indicating failure."""
    S = ar.compute_orientational_order_parameter(frags, verbose=False)
    return S < 0.6


def check_hop(frags_latest, frags_prev, sim_dir):
    """Returns a list of fragments that meet hop criteria for spinoff."""
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

    spin_offs.sort(key=lambda x: -x[1])  # Sort by distance descending
    return spin_offs[:MAX_SPINOFFS_PER_STEP]


def write_focus_file(sim_dir, frag_index, from_col, dist, frame_number=None, to_col=None):
    """Writes molecule.focus with hopping info for a tracked fragment."""
    focus_path = sim_dir / "molecule.focus"
    with open(focus_path, "w") as f:
        f.write("Tracked Fragment Hopping Info\n")
        f.write("-----------------------------\n")
        f.write(f"Fragment index:     {frag_index}\n")
        f.write(f"From column:        {from_col}\n")
        if to_col is not None:
            f.write(f"To column:          {to_col}\n")
        f.write(f"Trigger distance:   {dist:.4f} A\n")
        if frame_number is not None:
            f.write(f"Detected at frame:  {frame_number}\n")


def check_focus(sim_dir, frags):
    """Checks if a tracked fragment has completed a hop or failed to do so."""
    focus_file = sim_dir / "molecule.focus"
    if not focus_file.exists():
        return None

    with open(focus_file) as f:
        lines = f.readlines()
        frag_index = int(lines[2].split(":")[1])
        orig_col = int(lines[3].split(":")[1])
        orig_dist = float(lines[5].split(":")[1].split()[0])  # Remove A

    coms = np.array([frag.center_of_mass() for frag in frags])
    box = frags[0].dimensions[:3]

    anchors, _ = ar.guess_column_positions(coms, box, return_noise=False)
    assigned_columns, r = ar.assign_fragments_to_columns(coms, anchors, box)
    current_col = assigned_columns[frag_index]

    anchor = anchors[orig_col]
    box_xy = box[:2]
    dist = np.linalg.norm((coms[frag_index][:2] - anchor + box_xy / 2) % box_xy - box_xy / 2)

    if current_col != orig_col:
        log_action(f"Success: Fragment {frag_index} hopped from {orig_col} to {current_col}")
        return 3
    elif dist > orig_dist:
        log_action(f"Re-spinoff: Fragment {frag_index} now {dist:.2f} Å from column {orig_col}")
        return 1
    elif dist < r * HOP_DISTANCE_FACTOR:
        log_action(f"Failed hop: Fragment {frag_index} returned to home column {orig_col}")
        return 2

    return None


def log_action(text):
    """Appends action text to log.txt with timestamp."""
    with open("log.txt", "a") as log:
        log.write(f"[{current_timestamp()}] {text}\n")


def current_timestamp():
    return str(np.datetime64('now')).replace("T", " ")
