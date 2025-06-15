import random

def analysis(xyz_folder_path):
    r = random.random()
    if r < 0.05:
        print("[ANALYSIS] Returning 2 (stop failed)")
        return 2
    elif r < 0.10:
        print("[ANALYSIS] Returning 3 (stop success)")
        return 3
    elif r < 0.30:
        print("[ANALYSIS] Returning 1 (spin off)")
        return 1
    else:
        print("[ANALYSIS] Returning 4 (continue)")
        return 4
