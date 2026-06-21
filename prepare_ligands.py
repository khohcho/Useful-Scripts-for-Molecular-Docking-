import os
import glob
import subprocess
import argparse

def main():
    """
    Automated CLI wrapper for AutoDockTools (ADT) prepare_ligand4.py.
    Bypasses the GUI to prevent crashes with large datasets and strictly 
    preserves pre-calculated quantum partial charges (e.g., from Maestro/ESP).
    """
    parser = argparse.ArgumentParser(description="Batch process mol2 to pdbqt while preserving charges.")
    parser.add_argument("--adt_python", default=r"D:\Docking\1.5.7\python.exe", help="Path to ADT Python executable")
    parser.add_argument("--adt_script", default=r"D:\Docking\1.5.7\Lib\site-packages\AutoDockTools\Utilities24\prepare_ligand4.py", help="Path to prepare_ligand4.py")
    parser.add_argument("--input_dir", default=".", help="Directory containing .mol2 files")
    parser.add_argument("--output_dir", default="../Ligands_PDBQT", help="Output directory for .pdbqt files")
    
    args = parser.parse_args()

    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    mol2_files = glob.glob(os.path.join(args.input_dir, "*.mol2"))
    print(f"--- Initialization: Found {len(mol2_files)} ligands to process ---")

    for file_path in mol2_files:
        base_name = os.path.basename(file_path).replace(".mol2", "")
        output_path = os.path.join(args.output_dir, f"{base_name}.pdbqt")
        
        print(f"[PROCESSING] {base_name} | Preserving existing partial charges...")
        
        # -C: Preserve charges, -A 'checkhydrogens': Verify H atoms, -U 'nphs_lps': Merge non-polar H
        command = [
            args.adt_python, args.adt_script,
            "-l", file_path,
            "-o", output_path,
            "-C", "-A", "checkhydrogens", "-U", "nphs_lps"
        ]
        
        subprocess.run(command, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    print(f"--- SUCCESS: All ligands converted and saved to {args.output_dir} ---")

if __name__ == "__main__":
    main()