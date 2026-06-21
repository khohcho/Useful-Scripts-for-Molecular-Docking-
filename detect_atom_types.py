import os
import glob
import argparse

def main():
    """
    Scans a directory of PDBQT files to identify all unique atom types.
    Critical for generating accurate AutoDock-GPU .dpf template files, 
    preventing runtime 'unrecognized atom type' fatal errors.
    """
    parser = argparse.ArgumentParser(description="Extract unique atom types from PDBQT library.")
    parser.add_argument("--input_dir", default=".", help="Directory containing ligand .pdbqt files")
    args = parser.parse_args()

    print(f"--- Scanning PDBQT files in '{args.input_dir}' ---")

    pdbqt_files = glob.glob(os.path.join(args.input_dir, "*.pdbqt"))
    unique_atoms = set()

    if not pdbqt_files:
        print("ERROR: No .pdbqt files found in the specified directory.")
        return

    for filepath in pdbqt_files:
        try:
            with open(filepath, 'r') as f:
                for line in f:
                    if line.startswith("ATOM") or line.startswith("HETATM"):
                        parts = line.split()
                        if len(parts) > 2:
                            atom_type = parts[-1].strip()
                            unique_atoms.add(atom_type)
        except Exception as e:
            print(f"Read error at {filepath} -> {e}")

    print("\n" + "="*50)
    print("UNIQUE ATOM TYPES EXTRACTED (Ready for .dpf generation):")
    print("="*50)
    print(" ".join(sorted(unique_atoms)))
    print("="*50)
    print(f"Total files parsed: {len(pdbqt_files)}")

if __name__ == "__main__":
    main()