import os
import glob
import math
import shutil

# --- ADVANCED CONFIGURATION ---
OUTPUT_DIR = "./output"
LOG_DIR = "./logs"
ANCHOR_ATOM_TYPE = "P"

# Primary Binding Site (e.g., Catalytic center / Cofactors)
PRIMARY_TARGET_COORDS = [(97.212, 88.488, 100.157), (92.442, 84.162, 100.658)]
THRESHOLDS = [4.0, 5.0, 6.0] 

# Secondary Binding Site (e.g., Allosteric pocket or secondary pharmacophore requirement)
SECONDARY_TARGET_COORD = (88.218, 96.291, 109.459)
SECONDARY_THRESH = 6.0 

# Filtering Modes
ONLY_UNIQUE_LIGANDS = False
TOP_N_LIMIT = 50

def calculate_distance(p1, p2):
    return math.sqrt(sum((a - b) ** 2 for a, b in zip(p1, p2)))

def analyze_multitarget_pose(pdbqt_file, threshold):
    fname = os.path.basename(pdbqt_file)
    
    # Read affinity from Vina Log
    log_name = f"{fname.split('_out')[0]}.log".replace("__", "_")
    log_path = os.path.join(LOG_DIR, log_name)
    energy = None
    if os.path.exists(log_path):
        with open(log_path, 'r') as f:
            for line in f:
                parts = line.split()
                if len(parts) >= 2 and parts[0] == "1":
                    energy = float(parts[1])
                    break
    if energy is None: return False, None, None

    all_coords = []
    anchor_coords = []
    try:
        with open(pdbqt_file, 'r') as f:
            for line in f:
                if "ENDMDL" in line: break
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    coords = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
                    all_coords.append(coords)
                    if line[77:79].strip() == ANCHOR_ATOM_TYPE:
                        anchor_coords.append({"id": int(line[6:11].strip()), "coords": coords})
    except: return False, None, None

    # 1. Secondary Target Verification
    target_passed = False
    min_secondary_dist = 999.0
    for c in all_coords:
        d = calculate_distance(c, SECONDARY_TARGET_COORD)
        if d < min_secondary_dist: min_secondary_dist = d
        if d <= SECONDARY_THRESH: target_passed = True
            
    if not target_passed or not anchor_coords:
        return False, None, None 

    # 2. Primary Target Verification
    anchor_coords.sort(key=lambda x: x["id"])
    primary_passed = False
    dist_details = []
    
    for idx, p in enumerate(anchor_coords[:3]):
        dists = [calculate_distance(p["coords"], t) for t in PRIMARY_TARGET_COORDS]
        if dists[0] <= threshold or dists[1] <= threshold:
            primary_passed = True
        dist_details.append(f"Anchor-{idx+1} -> Primary_1: {dists[0]:.3f}A, Primary_2: {dists[1]:.3f}A")

    if primary_passed:
        dist_details.append(f"Secondary Target -> Closest atom distance: {min_secondary_dist:.3f}A (PASSED)")
        return True, energy, "\n".join(dist_details)
    return False, None, None

def run_multitarget_analysis():
    pdbqt_files = glob.glob(os.path.join(OUTPUT_DIR, "*.pdbqt"))
    mode_info = f"TOP {TOP_N_LIMIT} UNIQUE" if ONLY_UNIQUE_LIGANDS else "ALL VALID POSES"
    print(f"--- Analyzing {len(pdbqt_files)} files [V3: Multi-Target Filter | Mode: {mode_info}] ---")

    for thr in THRESHOLDS:
        folder_name = f"Analysis_MultiTarget_{int(thr)}A"
        os.makedirs(folder_name, exist_ok=True)
        
        results = []
        for f in pdbqt_files:
            is_ok, energy, dist_text = analyze_multitarget_pose(f, thr)
            if is_ok:
                results.append({"filename": os.path.basename(f), "energy": energy, "dist_text": dist_text})

        results.sort(key=lambda x: x["energy"])
        report_path = os.path.join(folder_name, f"MultiTarget_Report_{int(thr)}A.txt")
        
        with open(report_path, "w", encoding="utf-8") as report:
            report.write(f"ADVANCED DOCKING REPORT - {mode_info} - THRESHOLD: {thr}A\n" + "="*85 + "\n\n")
            
            written_ligands = set()
            count, idx = 0, 0
            
            while idx < len(results):
                if ONLY_UNIQUE_LIGANDS and count >= TOP_N_LIMIT: break
                    
                res = results[idx]
                base_name = res["filename"].split('_out')[0]
                
                if ONLY_UNIQUE_LIGANDS and base_name in written_ligands:
                    idx += 1
                    continue
                
                report.write(f"{count+1}. LIGAND: {res['filename']}\n")
                report.write(f"Affinity: {res['energy']:.3f} kcal/mol\n")
                report.write(f"Spatial Constraints Verification:\n{res['dist_text']}\n")
                report.write("-" * 75 + "\n\n")
                
                count += 1
                written_ligands.add(base_name)
                idx += 1

def extract_multitarget_ligands():
    print("\nExtracting successful multi-target poses for downstream analysis...")
    for thr in THRESHOLDS:
        thr_int = int(thr)
        report_path = os.path.join(f"Analysis_MultiTarget_{thr_int}A", f"MultiTarget_Report_{thr_int}A.txt")
        target_dir = os.path.join(OUTPUT_DIR, f"Extracted_MultiTarget_{thr_int}A")
        
        if os.path.exists(report_path):
            os.makedirs(target_dir, exist_ok=True)
            copied_count = 0
            with open(report_path, "r", encoding="utf-8") as f:
                for line in f:
                    if ". LIGAND: " in line:
                        fname = line.split(". LIGAND: ")[1].strip()
                        src = os.path.join(OUTPUT_DIR, fname)
                        dst = os.path.join(target_dir, fname)
                        if os.path.exists(src):
                            shutil.copy2(src, dst)
                            copied_count += 1
            print(f" -> {copied_count} ligand(s) successfully copied to {target_dir}")

if __name__ == "__main__":
    run_multitarget_analysis()
    extract_multitarget_ligands()