from hf2.coordinator import SimulationManager

sim_root = "your/path/here"
manager = SimulationManager(sim_root, verbose = True)
manager.run(interval = 30)