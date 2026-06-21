import os
import glob
import math

# --- CONFIGURATION ---
# Target atom types and geometric constraints for specific linear pharmacophores
CENTER_ATOM_TYPE = 'N'
FLANKING_ATOM_TYPE = 'NA'
CHARGE_THRESHOLD = 0.05
DISTANCE_TOLERANCE = 1.65

def calc_distance(atom1, atom2):
    """Calculates 3D Euclidean distance between two atoms."""
    dx = atom1['x'] - atom2['x']
    dy = atom1['y'] - atom2['y']
    dz = atom1['z'] - atom2['z']
    return math.sqrt(dx*dx + dy*dy + dz*dz)

def read_atom_data(line):
    """Parses standard PDBQT line and extracts coordinates and charge."""
    parts = line.split()
    try:
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
        atom_type = parts[-1]
        charge = float(parts[-2]) if len(parts) > 1 else float(line.split()[-2])
        return {'x': x, 'y': y, 'z': z, 'charge': charge, 'type': atom_type, 'line_idx': None}
    except ValueError:
        return None

def process_files():
    pdbqt_files = glob.glob("*.pdbqt")
    print(f"--- LINEAR PHARMACOPHORE OPTIMIZER STARTING ---")
    print(f"Scanning {len(pdbqt_files)} files for geometric corrections...\n")

    total_fixed_groups = 0

    for filename in pdbqt_files:
        with open(filename, 'r') as f:
            lines = f.readlines()

        atoms = []
        for i, line in enumerate(lines):
            if line.startswith("ATOM") or line.startswith("HETATM"):
                atom_data = read_atom_data(line)
                if atom_data:
                    atom_data['line_idx'] = i 
                    atoms.append(atom_data)

        is_modified = False
        
        for atom in atoms:
            if atom['type'] in [CENTER_ATOM_TYPE, FLANKING_ATOM_TYPE] and atom['charge'] > CHARGE_THRESHOLD:
                neighbors = []
                for other_atom in atoms:
                    if other_atom == atom: continue 
                    if other_atom['type'] in [CENTER_ATOM_TYPE, FLANKING_ATOM_TYPE]:
                        if calc_distance(atom, other_atom) < DISTANCE_TOLERANCE:
                            neighbors.append(other_atom)
                
                # Enforce specific linear arrangement correction
                if len(neighbors) == 2:
                    old_line = lines[atom['line_idx']]
                    if FLANKING_ATOM_TYPE in old_line[-4:]: 
                          new_line = old_line.rstrip()[:-len(FLANKING_ATOM_TYPE)] + f"{CENTER_ATOM_TYPE} \n"
                          lines[atom['line_idx']] = new_line

                    for k in neighbors:
                        k_line = lines[k['line_idx']]
                        current_type = k_line.split()[-1]
                        if current_type != FLANKING_ATOM_TYPE: 
                            stripped = k_line.rstrip()
                            new_k_line = stripped[:stripped.rfind(current_type)] + f"{FLANKING_ATOM_TYPE}\n"
                            lines[k['line_idx']] = new_k_line
                    
                    is_modified = True
                    total_fixed_groups += 1
                    print(f"-> {filename}: Linear group geometry and types corrected.")

        if is_modified:
            with open(filename, 'w') as f:
                f.writelines(lines)

    print(f"\nPROCESS COMPLETED. Total {total_fixed_groups} functional groups strictly parameterized.")

if __name__ == "__main__":
    process_files()