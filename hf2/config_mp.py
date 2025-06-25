
# --------------------------------------------------
# Directory and File Naming
# --------------------------------------------------

# Prefix used to identify the reference .xyz file (e.g., "F7_ramp_original.xyz")
REF_PREFIX = "F8_ramp"

# Subdirectory inside each simulation path where collected files are stored
XYZ_SUBDIR = "yield"



# --------------------------------------------------
# Misc Control Settings
# --------------------------------------------------

# Normalize atom names (e.g., CA -> C, HA -> H) to be compatible with MDAnalysis
NORMALIZE_ATOM_NAMES = True

TINKER_INACTIVITY_TIMEOUT = 3600  # seconds (60 minutes)

# --------------------------------------------------
# TINKER information
# --------------------------------------------------

# Path to the compiled TINKER binary (can be relative to home or absolute)
TINKER_BINARY = "~/Desktop/TinkerExp/bin/"

TINKER_PARAMS = "~/Desktop/TinkerExp/params"

# --------------------------------------------------
# Monitor settings
# --------------------------------------------------

# Minimum number of output frames between checks
MIN_OBSERVATION_FRAMES = 5

# Minimum amount of wall time between checks (seconds)
MIN_OBSERVATION_TIME = 300

