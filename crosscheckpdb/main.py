# pip install biopython --upgrade  # https://biopython.org/wiki/Download
import argparse
from Bio import PDB
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument(
    '--path_before',
    type=str,
    default="/home/qiliu02/GHDDI/GHDDI-Database/dev_query_pdb/comp_before.pdb"
)
parser.add_argument(
    '--path_after',
    type=str,
    default="/home/qiliu02/GHDDI/GHDDI-Database/dev_query_pdb/comp.pdb"
)
parser.add_argument(
    '--query_str',
    type=str,
    default="Atom_Name == 'OD1' and Amino_Acid == 'ASP' and Chain_ID=='R' and Residue_Sequence_Number == '121'"
)
args = parser.parse_args()

def pdb_to_dataframe(pdb_file_path):
    # Create a PDB parser
    parser = PDB.PDBParser(QUIET=True)

    # Parse the PDB file
    structure = parser.get_structure("ATOM", pdb_file_path)

    # Create lists to store data
    atom_data = []

    # Iterate through the structure and collect atom coordinates
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    atom_info = {
                        "Record_Type": atom.get_full_id()[0],
                        "Atom_Serial_Number": atom.get_serial_number(),
                        "Atom_Name": atom.get_name(),
                        "Amino_Acid": residue.resname,
                        "Chain_ID": chain.id,
                        "Residue_Sequence_Number": residue.id[1],
                        "X_Coordinate": atom.coord[0],
                        "Y_Coordinate": atom.coord[1],
                        "Z_Coordinate": atom.coord[2],
                        "Occupancy": atom.get_occupancy(),
                        "B_Factor": atom.get_bfactor(),
                        "Element": atom.element,
                    }
                    atom_data.append(atom_info)

    # Create a DataFrame from the collected data
    df = pd.DataFrame(atom_data)

    return df

path_before = args.path_before # "/home/qiliu02/GHDDI/GHDDI-Database/dev_query_pdb/comp_before.pdb"
df_before = pdb_to_dataframe(path_before)
print(f"df_before: {df_before.shape}")

path_after = args.path_after # "/home/qiliu02/GHDDI/GHDDI-Database/dev_query_pdb/comp.pdb"
df_after = pdb_to_dataframe(path_after)
print(f"df_after: {df_after.shape}")

query_str = args.query_str # "Atom_Name == 'OD1' and Amino_Acid == 'ASP' and Chain_ID=='R' and Residue_Sequence_Number == '121'"
query_res = df_before.query(query_str)

try:
    print("... find in df_before:")
    print(query_res)
    x, y, z = query_res['X_Coordinate'].values[0], query_res['Y_Coordinate'].values[0], query_res['Z_Coordinate'].values[0]
    query_xyz = f"X_Coordinate == {x} and Y_Coordinate=={y} and Z_Coordinate=={z}"
    query_res_fin = df_after.query(query_xyz)
    print("... find in df_after:")
    print(query_res_fin)
except Exception as e:
    print(f"Error! {e}")
print("Finished.")
