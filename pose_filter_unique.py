import os
import glob
import math
import shutil

# --- CORE CONFIGURATION ---
OUTPUT_DIR = "./output"
LOG_DIR = "./logs"
PRIMARY_TARGET_COORDS = [(97.212, 88.488, 100.157), (92.442, 84.162, 100.658)]
ANCHOR_ATOM_TYPE = "P"
THRESHOLDS = [4.0, 5.0, 6.0] 
TOP_N_UNIQUE = 50

def calculate_distance(p1, p2):
    return math.sqrt(sum((a - b) ** 2 for a, b in zip(p1, p2)))

def get_affinity_data(pdbqt_filename):
    name_no_ext = pdbqt_filename.replace(".pdbqt", "")
    prefix = name_no_ext.split("_out")[0] if "_out" in name_no_ext else name_no_ext
    log_name = f"{prefix}.log".replace("__", "_")
    log_path = os.path.join(LOG_DIR, log_name)
    
    if not os.path.exists(log_path):
        return None
    try:
        with open(log_path, 'r', encoding="utf-8") as f:
            for line in f:
                parts = line.split()
                if len(parts) >= 2 and parts[0] == "1":
                    return float(parts[1])
    except: return None

def analyze_strict_m1(pdbqt_file, threshold):
    fname = os.path.basename(pdbqt_file)
    energy = get_affinity_data(fname)
    if energy is None: return False, None, None

    anchor_coords = []
    try:
        with open(pdbqt_file, 'r') as f:
            for line in f:
                if "ENDMDL" in line: break
                if (line.startswith("ATOM") or line.startswith("HETATM")) and line[77:79].strip() == ANCHOR_ATOM_TYPE:
                    coords = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
                    at_id = int(line[6:11].strip())
                    anchor_coords.append({"id": at_id, "coords": coords})
    except: return False, None, None

    if not anchor_coords: return False, None, None
    anchor_coords.sort(key=lambda x: x["id"])
    
    passed = False
    dist_details = []
    labels = ["Anchor-1", "Anchor-2", "Anchor-3"]
    
    for idx, p in enumerate(anchor_data[:3] if 'anchor_data' in locals() else anchor_coords[:3]):
        label = labels[idx] if idx < len(labels) else f"Anchor-{idx+1}"
        dists = [calculate_distance(p["coords"], target) for target in PRIMARY_TARGET_COORDS]
        if dists[0] <= threshold or dists[1] <= threshold:
            passed = True
        dist_details.append(f"{label} -> Target1: {dists[0]:.3f}A, Target2: {dists[1]:.3f}A")

    if passed:
        return True, energy, "\n".join(dist_details)
    return False, None, None

def run_analysis():
    pdbqt_files = glob.glob(os.path.join(OUTPUT_DIR, "*.pdbqt"))
    print(f"--- Processing {len(pdbqt_files)} files [V2: Unique Ligand Filter & Extraction] ---")

    for thr in THRESHOLDS:
        folder_name = f"Analysis_Unique_{int(thr)}A"
        os.makedirs(folder_name, exist_ok=True)
        
        results = []
        for f in pdbqt_files:
            is_ok, energy, dist_text = analyze_strict_m1(f, thr)
            if is_ok:
                results.append({"filename": os.path.basename(f), "energy": energy, "dist_text": dist_text})

        results.sort(key=lambda x: x["energy"])
        report_path = os.path.join(folder_name, f"Top_Unique_Report_{int(thr)}A.txt")
        
        with open(report_path, "w", encoding="utf-8") as report:
            report.write(f"DOCKING REPORT - TOP {TOP_N_UNIQUE} UNIQUE LIGANDS - THRESHOLD: {thr}A\n" + "="*85 + "\n\n")
            
            written_ligands = set()
            count, idx = 0, 0
            
            while count < TOP_N_UNIQUE and idx < len(results):
                res = results[idx]
                base_name = res["filename"].split('_out')[0]
                
                report.write(f"{idx+1}. LIGAND: {res['filename']}\n")
                report.write(f"Affinity: {res['energy']:.3f} kcal/mol\n")
                report.write(f"Proximity:\n{res['dist_text']}\n")
                report.write("-" * 75 + "\n\n")
                
                if base_name not in written_ligands:
                    count += 1
                    written_ligands.add(base_name)
                idx += 1

def extract_successful_ligands():
    print("\nExtracting successful binding poses for MD simulations...")
    for thr in THRESHOLDS:
        thr_int = int(thr)
        report_path = os.path.join(f"Analysis_Unique_{thr_int}A", f"Top_Unique_Report_{thr_int}A.txt")
        target_dir = os.path.join(OUTPUT_DIR, f"Extracted_{thr_int}A")
        
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
    run_analysis()
    extract_successful_ligands()