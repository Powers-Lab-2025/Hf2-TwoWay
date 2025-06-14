# Hf2

Hf2 is a flexible Python package for steering TINKER-based molecular dynamics simulations using Forward Flux Sampling (FFS). It continuously monitors simulation directories, converts `.dyn` files to `.xyz` format, and uses a user-defined analysis function to determine how to handle each simulation path.

## Overview

Each simulation path (e.g. A5-10-2) is managed by a `SimulationPath` object, which handles file conversion and calls an external `analysis()` function. Based on the result of the analysis, the simulation can be spun off, marked as successful, terminated, or allowed to continue. The top-level controller, `SimulationManager`, supervises all paths in a given directory, updating them periodically.

## Features

- Monitors `.dyn` files in active paths
- Converts `.dyn` files to `.xyz` format using a reference structure
- Supports flexible naming conventions and directory structures
- User-defined analysis function to guide path behavior
- Automatically spawns and initializes new simulation paths
- Marks stopped paths as either successful (V prefix) or failed (X prefix)
- Central manager tracks all paths and logs runtime activity

## Folder Structure

A typical simulation tree looks like:

```
SimSource/
├── A5/
│   ├── F7_ramp.001
│   ├── F7_ramp.002
│   ├── F7_ramp_original.xyz
│   └── XYZs/
├── A5-1/
│   └── ...
├── V5-2/          # Successful
├── X3-1/          # Failed
```

## Installation

Clone the repository and make sure the root directory is on your Python path:

```
git clone https://github.com/Powers-Lab-2025/Hf2.git
cd Hf2
```

## Usage

A basic run script can be found in `scripts/run_steering.py`. You can run it like:

```python
from hf2.coordinator import SimulationManager

manager = SimulationManager("/path/to/SimSource", verbose=True)
manager.run(interval=5)
```

This will initialize all active paths in the directory and update them every 5 seconds.

## Configuration

Behavior is controlled through `hf2/config.py`. You should define:

- `TINKER_START_COMMAND`: how to start a simulation
- `TINKER_STOP_COMMAND`: how to stop a simulation
- `REF_XYZ_PREFIX`: name of the reference .xyz file
- `XYZ_SUBDIR`: name of the directory to store converted xyz files
- `DYN_SUFFIX`: file suffix to match `.dyn` files

## Analysis

The package calls `hf2.analysis.analysis()` on each path’s most recent `.xyz` file. This function should return:

- `1` → spin off the path
- `2` → stop the path and mark as failed (X)
- `3` → stop the path and mark as success (V)
- `4` → continue running

The analysis function is entirely user-defined and should be edited to fit your system.

## License

MIT License
