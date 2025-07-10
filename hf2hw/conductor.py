import os
import time
import shutil
import numpy as np
from pathlib import Path
from datetime import datetime
from hf2.analysis import analysis
import hf2.config

class SimulationPath:
    def __init__(self, path, verbose=True):
        self.path = Path(path).resolve()
        self.verbose = verbose
        self.label = self.path.name
        self.frame_counter = 0
        self.last_activity_time = time.time()

        self.analysis_frame_counter = 0
        self.last_global_spin_frame = -hf2.config.COOLDOWN_GLOBAL
        self.molecule_last_spin_frame = {}

        if not (self.path / f"{hf2.config.REF_XYZ_PREFIX}.xyz").exists():
            raise FileNotFoundError(f"Missing reference xyz file: {hf2.config.REF_XYZ_PREFIX}.xyz")

        if self.verbose:
            print(f"[PASSIVE] Monitoring {self.label} in passive mode")

    def check_for_new_frames(self):
        """Check for new .xyz frames and update activity time"""
        xyz_files = sorted([
            f for f in self.path.iterdir()
            if f.is_file()
            and f.name.startswith(f"{hf2.config.REF_XYZ_PREFIX}.")
            and len(f.suffixes) == 1
            and f.suffix[1:].isdigit()
            and not f.name.endswith(".dyn")
            and not f.name.endswith(".key")
        ], key=lambda f: f.stat().st_mtime)
        
        if xyz_files:
            latest_mtime = xyz_files[-1].stat().st_mtime
            if not hasattr(self, 'last_frame_mtime') or latest_mtime > self.last_frame_mtime:
                self.last_frame_mtime = latest_mtime
                self.last_activity_time = time.time()
                self.frame_counter += 1
                if self.verbose:
                    print(f"[MONITOR] Detected new XYZ frame: {xyz_files[-1].name}")
                return True
        return False

    def is_inactive(self):
        """Check if TINKER has been inactive too long"""
        return time.time() - self.last_activity_time > hf2.config.TINKER_INACTIVITY_TIMEOUT

    def run_analysis(self):
        """
        Calls user-defined analysis() and dispatches the result.
        Logs the logline and handles filename-based spinoffs.
        """
        try:
            code, logline, filename = analysis(self.path, self)

            if logline:
                self._log_to_file(logline)

            if self.verbose:
                print(f"[ANALYSIS] Path {self.label} returned: code={code}, filename={filename}")

            return self.take_action(code, filename)

        except Exception as e:
            err_msg = f"[ERROR] Analysis failed for {self.label}: {e}"
            self._log_to_file(err_msg)
            if self.verbose:
                print(err_msg)
            return True

    def take_action(self, action_code, filename=""):
        if action_code == 1:
            self.spin_off(filename)
            return True
        elif action_code == 2:
            self.stop_as_failed()
            return False
        elif action_code == 3:
            self.stop_as_success()
            return False
        elif action_code == 4:
            self.continue_running()
            return True
        else:
            if self.verbose:
                print(f"[ACTION] Unknown or no-op action ({action_code}) for {self.label}")
            return True

    def spin_off(self, dyn_filename):
        """
        Creates a new simulation path directory using the given .dyn filename,
        the latest .xyz, and .key file.
        """
        base_name = self.label

        if hf2.config.PASSIVE_MODE and hf2.config.NESTED_SPINOFF_DIRS:
            spinoff_root = self.path / "spinoffs"
            spinoff_root.mkdir(exist_ok=True)
        else:
            spinoff_root = self.path.parent

        siblings = [d.name for d in spinoff_root.iterdir() if d.is_dir() and d.name.startswith(base_name + "-")]
        suffixes = [int(d.split("-")[-1]) for d in siblings if d.split("-")[-1].isdigit()]
        new_suffix = max(suffixes + [0]) + 1
        new_label = f"{base_name}-{new_suffix}"
        new_path = spinoff_root / new_label
        new_path.mkdir()

        ref_xyz_file = self.path / f"{hf2.config.REF_XYZ_PREFIX}.xyz"
        key_file = self.path / f"{hf2.config.REF_XYZ_PREFIX}.key"
        dyn_file = self.path / f"{hf2.config.REF_XYZ_PREFIX}.dyn"
        log_file = self.path / "log.txt"
        frame_xyz_file = self.path/f"{hf2.config.REF_XYZ_PREFIX}.{self.frame_counter}"

        to_copy = []

        if ref_xyz_file.exists():
            to_copy.append(ref_xyz_file)

        if dyn_file.exists():
            to_copy.append(dyn_file)

        if key_file.exists():
            to_copy.append(key_file)

        if log_file.exists():
            to_copy.append(log_file)

        if frame_xyz_file.exists():
            to_copy.append(frame_xyz_file)

        for f in to_copy:
            shutil.copy(f, new_path / f.name)

        if self.verbose:
            print(f"[SPINOFF] Created {new_label} with {len(to_copy)} files")
            print(f"[PASSIVE] {new_label} created in passive mode - manual TINKER startup required")

    def stop_as_failed(self):
        new_name = self.label.replace("A", "X", 1)
        self._rename_and_stop(new_name)
        self._log_to_file(f"[STOP] {self.label} marked as failed (X).")
        if self.verbose:
            print(f"[STOP] {self.label} marked as failed (X).")

        new_path = self.path.parent / new_name
        try:
            self.path.rename(new_path)
            self.path = new_path
            self.label = new_name
            if self.verbose:
                print(f"[RENAME] Directory renamed to {new_name}")
        except Exception as e:
            if self.verbose:
                print(f"[RENAME] Failed to rename directory: {e}")

    def stop_as_success(self):
        new_name = self.label.replace("A", "V", 1)
        self._rename_and_stop(new_name)
        self._log_to_file(f"[STOP] {self.label} marked as success (V).")
        if self.verbose:
            print(f"[STOP] {self.label} marked as success (V).")

        new_path = self.path.parent / new_name
        try:
            self.path.rename(new_path)
            self.path = new_path
            self.label = new_name
            if self.verbose:
                print(f"[RENAME] Directory renamed to {new_name}")
        except Exception as e:
            if self.verbose:
                print(f"[RENAME] Failed to rename directory: {e}")

    def _rename_and_stop(self, new_name):
        prefix = hf2.config.TINKER_PREFIX
        stop_file = self.path / f"{prefix}.end"
        stop_file.touch()

        if self.verbose:
            print(f"[TINKER STOP] Created stop file: {stop_file}")

    def continue_running(self):
        if self.verbose:
            print(f"[CONTINUE] {self.label} continuing without changes.")

    def update(self):
        if not self.path.exists():
            if self.verbose:
                print(f"[EXTERNAL] Path {self.label} no longer exists (externally removed/renamed)")
            return False

        current_parent = self.path.parent
        for pattern in ['X', 'V']:
            potential_new_name = self.label.replace('A', pattern, 1)
            potential_new_path = current_parent / potential_new_name
            if potential_new_path.exists() and not self.path.exists():
                if self.verbose:
                    print(f"[EXTERNAL] Detected external rename: {self.label} -> {potential_new_name}")
                return False

        current_name = self.path.name
        if 'X' in current_name or 'V' in current_name:
            if self.verbose:
                print(f"[STOP] Path {current_name} already marked as finished, stopping monitoring")
            return False

        has_new_frame = self.check_for_new_frames()
        
        if self.is_inactive():
            self._log_to_file(f"[TIMEOUT] No new frames for {hf2.config.TINKER_INACTIVITY_TIMEOUT}s")
            if self.verbose:
                print(f"[TIMEOUT] {self.label} inactive, stopping")
            return False
            
        if has_new_frame:
            return self.run_analysis()
        return True

    def _log_to_file(self, text):
        # Only log if the directory still exists
        if not self.path.exists():
            return
            
        log_path = self.path / "log.txt"
        timestamp = str(np.datetime64('now')).replace("T", " ")
        with open(log_path, "a") as log:
            log.write(f"[{timestamp}] {text}\n")