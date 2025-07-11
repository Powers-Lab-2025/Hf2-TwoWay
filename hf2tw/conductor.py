import os
import time
import shutil
import numpy as np
from pathlib import Path
from datetime import datetime
from hf2tw.analysis import analysis
import hf2tw.config

class SimulationPath:
    def __init__(self, path, verbose=True):
        self.path = Path(path).resolve()
        self.verbose = verbose
        self.label = self.path.name
        self.frame_counter = 0
        self.last_activity_time = time.time()

        self.analysis_frame_counter = 0
        
        self.direction = None
        if self.path.parent.name in ['fro', 'to']:
            self.direction = self.path.parent.name
            
        self.prefix = f"{hf2tw.config.REF_XYZ_PREFIX}_{self.direction}"
        
        self.molecule_displacement_frames = {}

        if not (self.path / f"{self.prefix}.xyz").exists():
            raise FileNotFoundError(f"Missing reference xyz file: {self.prefix}.xyz")

        if self.verbose:
            direction_str = f" [{self.direction}]" if self.direction else ""
            print(f"[PASSIVE] Monitoring {self.label}{direction_str} in passive mode")

    def check_for_new_frames(self):
        xyz_files = sorted([
            f for f in self.path.iterdir()
            if f.is_file()
            and f.name.startswith(f"{self.prefix}.")
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
        return time.time() - self.last_activity_time > hf2tw.config.TINKER_INACTIVITY_TIMEOUT

    def run_analysis(self):
        try:
            code, logline, _ = analysis(self.path, self)

            if logline:
                self._log_to_file(logline)

            if self.verbose:
                print(f"[ANALYSIS] Path {self.label} returned: code={code}")

            return self.take_action(code)

        except Exception as e:
            err_msg = f"[ERROR] Analysis failed for {self.label}: {e}"
            self._log_to_file(err_msg)
            if self.verbose:
                print(err_msg)
            return True

    def take_action(self, action_code):
        if action_code == 2:
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

    def stop_as_failed(self):
        new_name = self.label.replace("A", "X", 1)
        self._create_stop_file()
        self._log_to_file(f"[STOP] {self.label} marked as failed (X).")
        if self.verbose:
            direction_str = f" [{self.direction}]" if self.direction else ""
            print(f"[STOP] {self.label}{direction_str} marked as failed (X).")

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
        self._create_stop_file()
        self._log_to_file(f"[STOP] {self.label} marked as success (V).")
        if self.verbose:
            direction_str = f" [{self.direction}]" if self.direction else ""
            print(f"[STOP] {self.label}{direction_str} marked as success (V).")

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

    def _create_stop_file(self):
        stop_file = self.path / f"{self.prefix}.end"
        stop_file.touch()
        if self.verbose:
            print(f"[TINKER STOP] Created stop file: {stop_file}")

    def continue_running(self):
        if self.verbose:
            direction_str = f" [{self.direction}]" if self.direction else ""
            print(f"[CONTINUE] {self.label}{direction_str} continuing without changes.")

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
            self._log_to_file(f"[TIMEOUT] No new frames for {hf2tw.config.TINKER_INACTIVITY_TIMEOUT}s")
            if self.verbose:
                print(f"[TIMEOUT] {self.label} inactive, stopping")
            return False
            
        if has_new_frame:
            return self.run_analysis()
        return True

    def _log_to_file(self, text):
        if not self.path.exists():
            return
            
        log_path = self.path / "log.txt"
        timestamp = str(np.datetime64('now')).replace("T", " ")
        with open(log_path, "a") as log:
            log.write(f"[{timestamp}] {text}\n")