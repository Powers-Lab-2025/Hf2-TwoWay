# Hf2

Hf2 is a Python package for managing TINKER molecular dynamics simulations using Forward Flux Sampling. It monitors simulation directories, analyzes trajectory data, and automatically manages simulation branching based on user-defined criteria.

## Overview

The package works by monitoring active simulation paths and running analysis on new trajectory frames. Each simulation path is managed by a SimulationPath object that tracks frame generation and calls your analysis function. Based on the analysis results, simulations can continue, spawn new paths, or terminate. The SimulationManager coordinates multiple paths simultaneously.

## Features

- Monitors TINKER trajectory output files automatically
- User-defined analysis function determines simulation behavior
- Automatic spinoff creation when interesting events are detected
- Configurable cooldown periods to prevent excessive branching
- Inactivity detection to stop stalled simulations
- Supports both active and passive monitoring modes
- Flexible directory structure with nested spinoff organization

## Directory Structure

A typical simulation setup:

```
SimSource/
├── A1/
│   ├── F8_ramp.xyz
│   ├── F8_ramp.key
│   ├── F8_ramp.xyz
│   ├── log.txt
│   ├── F8_ramp.001
│   ├── F8_ramp.002
│   └── spinoffs/
│       ├── A1-1/
│       └── A1-2/
├── A2/
├── V3/          # Completed successfully
├── X4/          # Failed
```

Active paths start with "A", successful paths are renamed to "V", and failed paths become "X".

## Installation

Clone the repository and add it to your Python path:

```
git clone https://github.com/Powers-Lab-2025/Hf2.git
cd Hf2
```

## Quick Start

```python
from hf2.coordinator import SimulationManager

# Start monitoring simulations
manager = SimulationManager("/path/to/SimSource", verbose=True)
manager.run(interval=5)
```

## Configuration

Edit hf2/config.py to customize behavior:

Core settings:
- MAX_ACTIVE: Maximum simultaneous simulation paths
- REF_XYZ_PREFIX: Base name for trajectory files (e.g. "F8_ramp")
- TINKER_INACTIVITY_TIMEOUT: Seconds before inactive paths are stopped
- PASSIVE_MODE: If True, only monitors without starting new TINKER processes

TINKER settings (when not in passive mode):
- TINKER_BINARY: Path to TINKER executable
- TINKER_NUM_STEPS: Simulation steps per run
- TINKER_TIMESTEP_FS: Timestep in femtoseconds
- TINKER_TEMP: Simulation temperature

## Analysis Function

The package calls analysis(sim_dir, path_state) for each new frame. Your analysis function should examine the trajectory and return:

- 1: Create a spinoff simulation
- 2: Stop and mark as failed
- 3: Stop and mark as successful  
- 4: Continue running

The analysis function receives the simulation directory path and a state object for tracking per-path variables like frame counters and cooldowns.

Example analysis structure:

```python
def analysis(sim_dir, path_state):
    # Read trajectory data
    # Perform your analysis
    # Check for interesting events
    
    if spinoff_condition_met:
        return 1, "Spinoff reason", "filename.dyn"
    elif failure_condition:
        return 2, "Failure reason", ""
    elif success_condition:
        return 3, "Success reason", ""
    else:
        return 4, "", ""
```

## Running Simulations

For passive monitoring of existing TINKER jobs:
1. Set PASSIVE_MODE = True in config.py
2. Point the manager at your simulation directory
3. The package will monitor and analyze without starting new processes

For active simulation management:
1. Set PASSIVE_MODE = False
2. Configure TINKER paths and parameters
3. The package will start and manage TINKER processes automatically

## License

MIT License