import pandas as pd
from rdkit import Chem
from rdkit.Chem import RDKFingerprint
from .base_fingerprint import Fingerprint

class RDKitFingerprint(Fingerprint):
    def calculate(self, input_file, output_file):
        df = pd.read_csv(input_file)
        if 'smiles' not in df.columns:
            raise ValueError("Input file must contain a 'smiles' column")
        df = df.drop_duplicates(subset=['smiles'])
        smiles_list = df['smiles'].tolist()
        
        def get_fingerprint(smiles):
            mol = Chem.MolFromSmiles(smiles)
            if mol is not None:
                return list(RDKFingerprint(mol))
            else:
                return None
        
        fingerprints = [get_fingerprint(smiles) for smiles in smiles_list]
        fingerprints_df = pd.DataFrame(fingerprints)
        result_df = pd.concat([df.reset_index(drop=True), fingerprints_df.reset_index(drop=True)], axis=1)
        result_df.to_csv(output_file, index=False)
        return result_df