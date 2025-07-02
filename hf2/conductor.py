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

        if not (self.path / f"{hf2.config.REF_XYZ_PREFIX}.xyz").exists():
            raise FileNotFoundError(f"Missing reference xyz file: {hf2.config.REF_XYZ_PREFIX}.xyz")

        if not hf2.config.PASSIVE_MODE:
            cmd = f"{hf2.config.TINKER_BINARY} {hf2.config.TINKER_PREFIX} {hf2.config.TINKER_NUM_STEPS} " \
                  f"{hf2.config.TINKER_TIMESTEP_FS} {hf2.config.TINKER_SNAPSHOT_INTERVAL_PS} " \
                  f"{hf2.config.TINKER_CONTROL_FLAGS} {hf2.config.TINKER_TEMP} 1 > " \
                  f"{hf2.config.TINKER_PREFIX}.{hf2.config.TINKER_RECORD_INDEX}.out"
        
            if self.verbose:
                print(f"[TINKER START] Running in {self.path}:\n  {cmd}")
        
            os.system(f"cd {self.path} && {cmd}")
        else:
            if self.verbose:
                print(f"[PASSIVE] Skipping TINKER startup in {self.label}")

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
            code, logline, filename = analysis(self.path)

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

        # Identify files to copy
        xyz_files = sorted(self.path.glob(f"{hf2.config.REF_XYZ_PREFIX}.*"), key=lambda f: f.stat().st_mtime)
        key_file = self.path / f"{hf2.config.REF_XYZ_PREFIX}.key"
        dyn_file = self.path / dyn_filename if dyn_filename else None

        to_copy = []
        if xyz_files:
            to_copy.append(xyz_files[-1])
        if dyn_file and dyn_file.exists():
            to_copy.append(dyn_file)
        if key_file.exists():
            to_copy.append(key_file)

        for f in to_copy:
            shutil.copy(f, new_path / f.name)

        if self.verbose:
            print(f"[SPINOFF] Created {new_label} with {len(to_copy)} files")

        if hf2.config.PASSIVE_MODE:
            if self.verbose:
                print(f"[PASSIVE] {new_label} created but not started (passive mode ON)")
        else:
            if self.verbose:
                print(f"[SPINOFF] Launching new TINKER instance in {new_label}")
            os.system(hf2.config.TINKER_START_COMMAND.format(path=new_path))

    def stop_as_failed(self):
        new_name = self.label.replace("A", "X", 1)
        self._rename_and_stop(new_name)
        self._log_to_file(f"[STOP] {self.label} marked as failed (X).")
        if self.verbose:
            print(f"[STOP] {self.label} marked as failed (X).")

    def stop_as_success(self):
        new_name = self.label.replace("A", "V", 1)
        self._rename_and_stop(new_name)
        self._log_to_file(f"[STOP] {self.label} marked as success (V).")
        if self.verbose:
            print(f"[STOP] {self.label} marked as success (V).")

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
        # Check for new frames first
        has_new_frame = self.check_for_new_frames()
        
        # Check for inactivity
        if self.is_inactive():
            self._log_to_file(f"[TIMEOUT] No new frames for {hf2.config.TINKER_INACTIVITY_TIMEOUT}s")
            if self.verbose:
                print(f"[TIMEOUT] {self.label} inactive, stopping")
            return False  # Signal to remove this path
            
        # Only run analysis if there's a new frame
        if has_new_frame:
            return self.run_analysis()
        return True  # Continue running

    def _log_to_file(self, text):
        log_path = self.path / "log.txt"
        timestamp = str(np.datetime64('now')).replace("T", " ")
        with open(log_path, "a") as log:
            log.write(f"[{timestamp}] {text}\n")