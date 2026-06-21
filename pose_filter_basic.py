import os
import glob
import math

# --- CORE CONFIGURATION ---
OUTPUT_DIR = "./output"
LOG_DIR = "./logs"
PRIMARY_TARGET_COORDS = [(97.212, 88.488, 100.157), (92.442, 84.162, 100.658)] # Receptors/Cofactors
ANCHOR_ATOM_TYPE = "P" 
THRESHOLDS = [4.0, 5.0, 6.0] 
BASE_TOP_N = 50

def calculate_distance(p1, p2):
    return math.sqrt(sum((a - b) ** 2 for a, b in zip(p1, p2)))

def get_precise_m1_energy(pdbqt_filename):
    clean_name = pdbqt_filename.replace("_out", "").replace(".pdbqt", "")
    if clean_name.endswith("_"): clean_name = clean_name[:-1]
    
    log_path = os.path.join(LOG_DIR, f"{clean_name}.log")
    if not os.path.exists(log_path): return None

    try:
        with open(log_path, 'r', encoding="utf-8") as f:
            for line in f:
                parts = line.split()
                if len(parts) >= 2 and parts[0] == "1":
                    return float(parts[1])
    except: return None
    return None

def analyze_pose_distance(pdbqt_file, threshold):
    fname = os.path.basename(pdbqt_file)
    m1_energy = get_precise_m1_energy(fname)
    if m1_energy is None: return False, None, None, None

    with open(pdbqt_file, 'r') as f:
        content = f.read()
    
    models = content.split("ENDMDL")[:4]
    if not models or len(models[0].strip()) == 0: return False, None, None, None
    
    valid_pose = None
    best_dist_text = ""
    
    for i, model_content in enumerate(models):
        anchor_data = [] 
        for line in model_content.splitlines():
            if (line.startswith("ATOM") or line.startswith("HETATM")) and line[77:79].strip() == ANCHOR_ATOM_TYPE:
                try:
                    at_id = int(line[6:11].strip())
                    coords = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
                    anchor_data.append({"id": at_id, "coords": coords})
                except: continue
        
        if not anchor_data: continue
        anchor_data.sort(key=lambda x: x["id"])
        
        model_passed = False
        temp_dist_lines = []
        labels = ["Anchor-1", "Anchor-2", "Anchor-3"]
        
        for idx, p in enumerate(anchor_data[:3]):
            label = labels[idx] if idx < len(labels) else f"Anchor-{idx+1}"
            d1 = calculate_distance(p["coords"], PRIMARY_TARGET_COORDS[0])
            d2 = calculate_distance(p["coords"], PRIMARY_TARGET_COORDS[1])
            temp_dist_lines.append(f"{label} -> Target1: {d1:.3f}A, Target2: {d2:.3f}A")
            if d1 <= threshold or d2 <= threshold:
                model_passed = True
        
        if model_passed:
            valid_pose = i + 1
            best_dist_text = "\n".join(temp_dist_lines)
            break 
            
    return (True, m1_energy, valid_pose, best_dist_text) if valid_pose else (False, None, None, None)

def run():
    pdbqt_files = glob.glob(os.path.join(OUTPUT_DIR, "*.pdbqt"))
    print(f"--- Analyzing {len(pdbqt_files)} docking outputs [V1: Basic Distance Filter] ---")

    for thr in THRESHOLDS:
        folder_name = f"Analysis_Basic_{int(thr)}A"
        os.makedirs(folder_name, exist_ok=True)
        
        results = []
        for f in pdbqt_files:
            is_ok, m1_energy, vize_poz, dist_text = analyze_pose_distance(f, thr)
            if is_ok:
                results.append({"filename": os.path.basename(f), "energy": m1_energy, "pose": vize_poz, "dist_text": dist_text})

        results.sort(key=lambda x: x["energy"])
        
        report_path = os.path.join(folder_name, f"Top_Filtered_Report_{int(thr)}A.txt")
        with open(report_path, "w", encoding="utf-8") as report:
            report.write(f"VIRTUAL SCREENING ANALYSIS - THRESHOLD: {thr}A\n")
            report.write("="*85 + "\n\n")
            
            written_ligands = set()
            count = 0
            idx = 0
            while count < BASE_TOP_N and idx < len(results):
                res = results[idx]
                base_name = res["filename"].split('_out')[0]
                
                report.write(f"{idx+1}. LIGAND: {res['filename']}\n")
                if res["pose"] > 1:
                    report.write(f"-> Selected Pose: {res['pose']}\n")
                
                report.write(f"Binding Affinity: {res['energy']:.3f} kcal/mol\n")
                report.write(f"Geometric Proximity:\n{res['dist_text']}\n")
                report.write("-" * 75 + "\n\n")
                
                if base_name not in written_ligands:
                    count += 1
                    written_ligands.add(base_name)
                idx += 1
        
        print(f"Completed processing for {thr}A threshold.")

if __name__ == "__main__":
    run()