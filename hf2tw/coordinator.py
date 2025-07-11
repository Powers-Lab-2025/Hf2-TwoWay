import time
import subprocess
import shutil
from pathlib import Path
from collections import deque
from datetime import datetime
from hf2tw.conductor import SimulationPath
import hf2tw.config

class SimulationManager:
    def __init__(self, sim_root, verbose=True):
        self.sim_root = Path(sim_root).resolve()
        self.paths = {}
        self.verbose = verbose
        self.start_time = datetime.now()

        if self.verbose:
            print(f"[MANAGER] Monitoring simulation root: {self.sim_root}")
        
        self.scan_for_new_paths()

    def scan_for_new_paths(self):
        for direction in ['fro', 'to']:
            direction_path = self.sim_root / direction
            
            if not direction_path.exists() or not direction_path.is_dir():
                if self.verbose:
                    print(f"[MANAGER] Directory {direction}/ not found, skipping...")
                continue
            
            for subdir in direction_path.iterdir():
                if not subdir.is_dir():
                    continue
                    
                label = subdir.name
                if not label.startswith("A"):
                    continue
                    
                path_key = f"{direction}/{label}"
                
                if path_key in self.paths:
                    continue
        
                try:
                    conductor = SimulationPath(subdir, verbose=self.verbose)
                    self.paths[path_key] = conductor
                    if self.verbose:
                        print(f"[MANAGER] Activated new path: {path_key}")
                except Exception as e:
                    print(f"[MANAGER ERROR] Could not initialize {path_key}: {e}")

    def update_all(self):
        to_remove = []
    
        for path_key, conductor in list(self.paths.items()):
            try:
                still_active = conductor.update()
                if not still_active:
                    to_remove.append(path_key)
            except Exception as e:
                print(f"[UPDATE ERROR] Failed to update {path_key}: {e}")
    
        for path_key in to_remove:
            if self.verbose:
                print(f"[MANAGER] Removing stopped path: {path_key}")
            del self.paths[path_key]

    def count_active_by_direction(self):
        fro_count = sum(1 for key in self.paths if key.startswith("fro/"))
        to_count = sum(1 for key in self.paths if key.startswith("to/"))
        return fro_count, to_count

    def get_next_label(self, direction):
        direction_path = self.sim_root / direction
        if not direction_path.exists():
            direction_path.mkdir(parents=True)
            return "A1"
        
        existing_dirs = [d.name for d in direction_path.iterdir() if d.is_dir()]
        existing_numbers = []
        for name in existing_dirs:
            if name.startswith("A") and name[1:].isdigit():
                existing_numbers.append(int(name[1:]))
            elif name.startswith("X") and name[1:].isdigit():
                existing_numbers.append(int(name[1:]))
            elif name.startswith("V") and name[1:].isdigit():
                existing_numbers.append(int(name[1:]))
        
        next_number = max(existing_numbers + [0]) + 1
        return f"A{next_number}"

    def launch_new_simulation(self, direction):
        new_label = self.get_next_label(direction)
        new_path = self.sim_root / direction / new_label
        new_path.mkdir(parents=True)
        
        template_dir = Path(hf2tw.config.TEMPLATE_DIR)
        prefix = f"{hf2tw.config.REF_XYZ_PREFIX}_{direction}"
        
        files_to_copy = [
            template_dir / f"{prefix}.xyz",
            template_dir / f"{prefix}.key",
            template_dir / f"{prefix}.dyn"
        ]
        
        for src_file in files_to_copy:
            if src_file.exists():
                shutil.copy(src_file, new_path / src_file.name)
            else:
                print(f"[ERROR] Template file not found: {src_file}")
                return False
        
        tinker_command = hf2tw.config.TINKER_COMMAND.format(prefix=prefix)
        
        try:
            process = subprocess.Popen(
                tinker_command.split(),
                cwd=new_path,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL
            )
            
            if self.verbose:
                print(f"[LAUNCH] Started new simulation {direction}/{new_label} with PID {process.pid}")
                print(f"[LAUNCH] Command: {tinker_command}")
            
            return True
            
        except Exception as e:
            print(f"[LAUNCH ERROR] Failed to start TINKER: {e}")
            return False

    def maintain_active_simulations(self):
        fro_count, to_count = self.count_active_by_direction()
        min_active = hf2tw.config.MIN_ACTIVE_SIMULATIONS
        
        if fro_count < min_active:
            launches_needed = min_active - fro_count
            if self.verbose:
                print(f"[AUTO-LAUNCH] Forward simulations: {fro_count}/{min_active}, launching {launches_needed} new")
            for _ in range(launches_needed):
                self.launch_new_simulation('fro')
        
        if to_count < min_active:
            launches_needed = min_active - to_count
            if self.verbose:
                print(f"[AUTO-LAUNCH] Backward simulations: {to_count}/{min_active}, launching {launches_needed} new")
            for _ in range(launches_needed):
                self.launch_new_simulation('to')

    def log_status(self):
        elapsed = datetime.now() - self.start_time
        
        fro_count, to_count = self.count_active_by_direction()
        
        print("\n=== Simulation Status ===")
        print(f"Uptime: {elapsed}")
        print(f"Active paths: {len(self.paths)} total")
        print(f"  - Forward (fro): {fro_count}")
        print(f"  - Backward (to): {to_count}")
        print("=========================\n")

    def run(self, interval=5):
        if self.verbose:
            print(f"[MANAGER] Starting main loop with {interval}s interval.")
            print(f"[MANAGER] Monitoring for two-way shooting: fro/ and to/ directories")
            print(f"[MANAGER] Auto-launch enabled: maintaining {hf2tw.config.MIN_ACTIVE_SIMULATIONS} active simulations per direction")
            
        while True:
            self.scan_for_new_paths()
            self.update_all()
            self.maintain_active_simulations()
            self.log_status()
            time.sleep(interval)