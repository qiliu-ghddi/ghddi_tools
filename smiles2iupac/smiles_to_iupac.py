import argparse
import sys
import time
import pandas as pd
import multiprocessing as mp
from datetime import datetime
from pathlib import Path
# from dockstring import load_target
from joblib import Parallel, delayed
from tqdm import tqdm
from rdkit import Chem  # conda install -c conda-forge rdkit
from rdkit.Chem import AllChem
from rdkit import DataStructs
from STOUT import translate_forward, translate_reverse
import warnings
warnings.filterwarnings("ignore")



def parse_arguments():
    parser = argparse.ArgumentParser(description="Your program description here.")
    parser.add_argument('-f', '--dataset_filename', type=str, default='BBBP100.csv',
                        help='Filename for the dataset (default: "BBBP100.csv")')    
    parser.add_argument('--smiles_col', type=str, default='smiles',
                        help='Column name for smiles data (default: "smiles")')
    parser.add_argument('--index_col', type=str, default='num',
                        help='Column name for index (default: "num")')
    parser.add_argument('--is_parallel', action='store_true', help='If using the multiprocessing')
    return parser.parse_args()



def are_same_molecule(smiles1, smiles2):
    mol1 = Chem.MolFromSmiles(smiles1)
    mol2 = Chem.MolFromSmiles(smiles2)
    
    if mol1 is None or mol2 is None:
        return False
    
    canonical_smiles1 = Chem.MolToSmiles(mol1)
    canonical_smiles2 = Chem.MolToSmiles(mol2)
    
    if canonical_smiles1 == canonical_smiles2:
        return True
    else:
        return False



def calculate_tanimoto_similarity(smiles1, smiles2):
    mol1 = Chem.MolFromSmiles(smiles1)
    mol2 = Chem.MolFromSmiles(smiles2)
    
    if mol1 is None or mol2 is None:
        return None
    
    fp1 = AllChem.GetMorganFingerprint(mol1, 2)
    fp2 = AllChem.GetMorganFingerprint(mol2, 2)
    
    similarity = DataStructs.TanimotoSimilarity(fp1, fp2)
    return similarity



def smiles_to_iupac_on_parameter(index, param):
    sample = param['sample']
    index = sample.name
       
    smiles = sample[smiles_col]
    t_begin = time.time()
    

    IUPAC_name = translate_forward(smiles)
    print("IUPAC name of "+smiles+" is: "+IUPAC_name)
    smiles_reverse = translate_reverse(IUPAC_name)
    print("SMILES of "+IUPAC_name+" is: "+smiles_reverse)
        
    t_end = time.time()
    time_elasped = t_end - t_begin
    print(f"""Coverting SMILES {smiles} to IUPAC {IUPAC_name}. 
          And then IUPAC {IUPAC_name} to SMILES {smiles_reverse}.
          Total time elasped: {time_elasped} s""")
    
    result = {
        index_col: index,
        "SMILES": smiles,
        "IUPAC_name": IUPAC_name,
        "SMILES_reverse": smiles_reverse,
        "are_same_molecule": are_same_molecule(smiles, smiles_reverse),
        "tanimoto_similarity": calculate_tanimoto_similarity(smiles, smiles_reverse),
        'time_elasped': time_elasped
    }
    return result


def smiles_to_iupac_(index, sample):
    smiles = sample[smiles_col]
    t_begin = time.time()
    
    try:
        IUPAC_name = translate_forward(smiles)
        print("IUPAC name of "+smiles+" is: "+IUPAC_name)
        smiles_reverse = translate_reverse(IUPAC_name)
        print("SMILES of "+IUPAC_name+" is: "+smiles_reverse)
        t_end = time.time()
        time_elasped = t_end - t_begin
        print(f"""{index} ... Coverting SMILES {smiles} to IUPAC {IUPAC_name}. 
            And then IUPAC {IUPAC_name} to SMILES {smiles_reverse}.
            Total time elasped: {time_elasped} s""")
        
        result = {
            index_col: index,
            "SMILES": smiles,
            "IUPAC_name": IUPAC_name,
            "SMILES_reverse": smiles_reverse,
            "are_same_molecule": are_same_molecule(smiles, smiles_reverse),
            "tanimoto_similarity": calculate_tanimoto_similarity(smiles, smiles_reverse),
            'time_elasped': time_elasped
        }
        return result
    except Exception as e:
        print(f'Error!!! {index},{smiles},{e}')
    return None


def smiles_to_iupac(df):
    processed_list = []
    for i, sample in df.iterrows():
        res = smiles_to_iupac_(i, sample)
        if res:
            processed_list.append(res)
    return processed_list
        

def main(args):
    global smiles_col, index_col, dataset_filename, save_dir, is_parallel
    smiles_col = args.smiles_col
    index_col = args.index_col
    dataset_filename = args.dataset_filename
    save_dir = args.save_dir
    is_parallel = args.is_parallel
    Path(save_dir).mkdir(parents=True, exist_ok=True)
    
    # load samples
    df = pd.read_csv(
        dataset_filename, 
        sep=",",
        index_col=index_col
    )  
    print(f"#SMILES: {df.shape[0]}")

    print(f"Converting ...")
    if is_parallel:
        num_cores = mp.cpu_count()
        print(f"num of cpus: {num_cores}")        
        param_dict = {}
        for i, (index, line) in enumerate(df.iterrows()):
            # print(line)
            param_dict[i] = {'sample': line}        
        processed_list = Parallel(n_jobs=int(num_cores/2))(delayed(smiles_to_iupac_on_parameter)(i, param) for i, param in param_dict.items())
    else:
        processed_list = smiles_to_iupac(df)
    print(f"Converting finished ~")
    
    # save result
    print(f"Saving results ...")
    result_df = pd.DataFrame(processed_list)
    result_filename = Path(f"{save_dir}/{Path(dataset_filename).stem}-smiles_to_iupac.csv")
    result_df.to_csv(result_filename, index=False)

    print(f"Finished ~")
    

if __name__ == "__main__":
    t_stamp = datetime.now().strftime("%Y%m%d-%H%M%S")
    # log
    orig_stdout = sys.stdout
    f = open(f'{t_stamp}.txt', 'w')
    sys.stdout = f

    args = parse_arguments()
    # main
    main(args)
    
    sys.stdout = orig_stdout
    f.close()

