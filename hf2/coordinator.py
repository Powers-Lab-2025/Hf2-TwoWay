import time
from pathlib import Path
from collections import deque
from datetime import datetime
from hf2.conductor import SimulationPath

class SimulationManager:
    def __init__(self, sim_root, verbose=True):
        """
        Initializes the simulation manager.

        Parameters:
            sim_root (str or Path): Root directory where simulation paths are stored.
            verbose (bool): Whether to print log messages.
        """
        self.sim_root = Path(sim_root).resolve()
        self.paths = {}
        self.verbose = verbose
        self.start_time = datetime.now()
        self.recent_spinoffs = deque(maxlen=5)

        if self.verbose:
            print(f"[MANAGER] Monitoring simulation root: {self.sim_root}")
        
        self.scan_for_new_paths()

    def scan_for_new_paths(self):
        """
        Scans the root directory for new active simulation paths.

        Paths must:
        - Be a directory
        - Start with 'A'
        - Not already be registered

        If a valid path is found, it creates a SimulationPath object and registers it.
        """
        for subdir in self.sim_root.iterdir():
            if not subdir.is_dir():
                continue
            label = subdir.name
            if label.startswith("A") and label not in self.paths:
                try:
                    conductor = SimulationPath(subdir, verbose=self.verbose)
                    self.paths[label] = conductor
                    if self.verbose:
                        print(f"[MANAGER] Registered new active path: {label}")
                except Exception as e:
                    print(f"[MANAGER ERROR] Could not initialize {label}: {e}")

    def update_all(self):
        """
        Calls .update() on all registered simulation paths.

        This checks for new .dyn files, runs analysis, and triggers actions as needed.
        """
        for label, conductor in list(self.paths.items()):
            try:
                conductor.update()
            except Exception as e:
                print(f"[UPDATE ERROR] Failed to update {label}: {e}")

    def log_status(self):
        """
        Logs summary information:
        - Uptime
        - Number of active paths
        - Recent spinoffs (most recent 5)
        """
        elapsed = datetime.now() - self.start_time
        print("\n=== Simulation Status ===")
        print(f"Uptime: {elapsed}")
        print(f"Active paths: {len(self.paths)}")
        print("Recent spinoffs:")
        for item in self.recent_spinoffs:
            print(f" - {item}")
        print("=========================\n")

    def run(self, interval=5):
        """
        Starts the main loop.

        Every interval seconds, it:
        - Scans for new paths
        - Updates all paths
        - Logs status

        Parameters:
            interval (int): Number of seconds to wait between updates.
        """
        if self.verbose:
            print(f"[MANAGER] Starting main loop with {interval}s interval.")
        while True:
            self.scan_for_new_paths()
            self.update_all()
            self.log_status()
            time.sleep(interval)
