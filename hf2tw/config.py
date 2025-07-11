# =============================================================================
# HF2 Configuration File
# =============================================================================
# This file contains all configuration parameters for the HF2. 
# Modify these values to customize behavior for
# your specific molecular dynamics simulations.

# =============================================================================
# File and Directory Settings
# =============================================================================

# Base prefix for all simulation files (without extensions)
# This should match your TINKER simulation setup
# Example: if your files are "F8_ramp.xyz", "F8_ramp.key", etc., use "F8_ramp"
REF_XYZ_PREFIX = "F8_ramp"

# Whether to create spinoff directories inside parent directories
# True: creates A1/spinoffs/A1-1/
# False: creates A1-1/ alongside A1/
NESTED_SPINOFF_DIRS = True

# What files to copy when creating spinoffs
# "minimal": copy only .xyz, .key, .dyn, and log.txt files
# "all": copy all files (not currently implemented)
SPINOFF_COPY_MODE = "minimal"

# =============================================================================
# Monitoring and Timeout Settings  
# =============================================================================

# How long to wait (in seconds) before considering a simulation inactive
# If no new trajectory frames appear within this time, the path will be stopped
TINKER_INACTIVITY_TIMEOUT = 300  # 5 minutes

# =============================================================================
# Analysis and Detection Settings (User defined things)
# =============================================================================

# How often (in analysis frames) to check for orientational order failures
# Set to 0 to disable orientational checks
ORIENT_FAIL_INTERVAL = 10

# How often (in analysis frames) to check for molecular hopping events
# This controls the frequency of spinoff detection
HOP_CHECK_INTERVAL = 4

# Distance threshold for detecting molecular hops
# Molecules are considered "hopping" if they move farther than:
# HOP_DISTANCE_FACTOR x average_inter_column_distance
HOP_DISTANCE_FACTOR = 0.4

# Maximum number of spinoffs to create in a single analysis step
# Helps prevent computational overload from too many simultaneous branches
MAX_SPINOFFS_PER_STEP = 2

# How many frames back to compare against for hop detection
# Instead of comparing to the immediately previous frame, compare to N frames ago
HOP_REFERENCE_FRAMES_BACK = 4

PASSIVE_MODE = True

# =============================================================================
# Cooldown Period Settings
# =============================================================================
# These prevent excessive spinoff creation by enforcing waiting periods

# Global cooldown: minimum frames to wait after ANY spinoff before creating another
# Applies to the entire simulation path
COOLDOWN_GLOBAL = 5

# Per-molecule cooldown: minimum frames to wait before the same molecule can spinoff again
# Prevents a single molecule from creating too many branches
COOLDOWN_MOLECULE = 20

# =============================================================================
# Auto-launch Settings
# =============================================================================

# Minimum number of active simulations to maintain per direction
MIN_ACTIVE_SIMULATIONS = 4

# TINKER command to execute when launching new simulations
# Use {prefix} as a placeholder for the actual prefix (e.g., F8_ramp_fro)
# Example: "dynamic {prefix} 1000000 0.5 1 1"
TINKER_COMMAND = "dynamic {prefix} 1000000 0.5 1 1"

# Template directory containing initial files for new simulations
# Should contain {prefix}_fro.xyz, .key, .dyn for forward direction
# and {prefix}_to.xyz, .key, .dyn for backward direction
TEMPLATE_DIR = "./templates"