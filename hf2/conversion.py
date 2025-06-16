import os
import re
import numpy as np
from pathlib import Path
import hf2.config

def _read_dyn_coordinates(dyn_path, verbose=False):
    with open(dyn_path, 'r') as file:
        lines = file.readlines()

    if len(lines) < 6:
        raise ValueError("Dyn file too short to contain data.")

    tmp = lines[1].split()
    num_atoms = int(tmp[0])
    title = lines[1][len(tmp[0])+1:] if len(tmp) > 1 else ""

    if verbose:
        print(f"[DYN] {dyn_path.name} — {num_atoms} atoms\nTitle: {title.strip()}")

    sdims = lines[3].split()
    dims = [float(s.split('D')[0]) * 10**float(s.split('D')[1]) for s in sdims]
    angs = [float(s.split('D')[0]) * 10**float(s.split('D')[1]) for s in lines[4].split()]
    coords = np.zeros((3, num_atoms), dtype=np.float64)
    for i in range(num_atoms):
        line = lines[6 + i].split()
        for j in range(3):
            base, exp = line[j].split('D')
            coords[j, i] = float(base) * 10**float(exp)

    return num_atoms, dims, angs, coords


def _get_reference_xyz(ref_path, num_atoms):
    with open(ref_path, 'r') as file:
        lines = file.readlines()

    ref_num_atoms = int(lines[0].strip().split()[0])
    if ref_num_atoms != num_atoms:
        raise ValueError(f"Atom count mismatch: dyn={num_atoms}, ref_xyz={ref_num_atoms}")

    header = lines[1]
    atom_lines = lines[2:2 + num_atoms]
    return header, atom_lines


def _write_xyz_file(out_path, num_atoms, header, atom_lines, coords):
    with open(out_path, 'w') as out:
        out.write(f"{num_atoms}\n{header}")
        for i in range(num_atoms):
            parts = atom_lines[i].strip().split()
            index = parts[0]
            atom_name = parts[1]

            if hf2.config.NORMALIZE_ATOM_NAMES:
                if atom_name == "HA":
                    atom_name = "H"
                elif atom_name == "CA":
                    atom_name = "C"

            mol_id = parts[5]
            atom_type = parts[6]
            bonds = parts[7:]
            bond_str = "  " + "  ".join(bonds) if bonds else ""
            out.write(f"{index:>5s}  {atom_name:<2s}  {coords[0,i]:12.6f}  {coords[1,i]:12.6f}  {coords[2,i]:12.6f}  {mol_id:>5s}  {atom_type:>5s}{bond_str}\n")


def convert_dyn_to_xyz(dyn_path, ref_prefix="F7_ramp_original", out_subdir="XYZs", verbose=False):
    """
    Convert a .dyn file into a bonded .xyz using a reference .xyz file for atom topology.

    Parameters
    ----------
    dyn_path : str or Path
        Path to the .dyn file.
    ref_prefix : str
        Prefix (not full path) of the reference .xyz file, expected to be in the same directory.
    out_subdir : str
        Subdirectory name to place output .xyz files (created if not exists).
    verbose : bool
        If True, prints detailed output.

    Returns
    -------
    out_path : Path
        Path to the written .xyz file.
    """
    dyn_path = Path(dyn_path)
    dyn_stem = dyn_path.stem
    base_dir = dyn_path.parent
    xyz_dir = base_dir / out_subdir
    xyz_dir.mkdir(exist_ok=True)

    num_atoms, dims, angs, coords = _read_dyn_coordinates(dyn_path, verbose)

    ref_xyz_path = base_dir / f"{ref_prefix}.xyz"
    if not ref_xyz_path.exists():
        raise FileNotFoundError(f"Reference .xyz not found: {ref_xyz_path}")

    header, atom_lines = _get_reference_xyz(ref_xyz_path, num_atoms)

    out_path = xyz_dir / f"{dyn_stem}_converted.xyz"
    _write_xyz_file(out_path, num_atoms, header, atom_lines, coords)

    if verbose:
        print(f"[XYZ] Wrote bonded .xyz to: {out_path}")

    return out_path
