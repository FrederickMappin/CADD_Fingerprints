import pandas as pd
from padelpy import from_smiles
from .base_fingerprint import Fingerprint

class PubchemFingerprint(Fingerprint):
    def calculate(self, input_file, output_file):
        # Read the input CSV file
        df = pd.read_csv(input_file)
        
        # Check if 'smiles' column exists
        if 'smiles' not in df.columns:
            raise ValueError("Input file must contain a 'smiles' column")
        
        # Remove duplicate SMILES strings
        df = df.drop_duplicates(subset=['smiles'])
        
        # Extract the list of SMILES strings
        smiles_list = df['smiles'].tolist()
        
        # Calculate fingerprints for the list of SMILES strings with threading
        fingerprints = from_smiles(smiles_list, fingerprints=True, descriptors=False, threads=4)
        
        # Convert the list of dictionaries to a DataFrame
        fingerprints_df = pd.DataFrame(fingerprints)
        
        # Concatenate the original DataFrame with the fingerprints DataFrame
        result_df = pd.concat([df.reset_index(drop=True), fingerprints_df.reset_index(drop=True)], axis=1)
        
        # Save the results to the output CSV file
        result_df.to_csv(output_file, index=False)
        
        # Return the resulting DataFrame
        return result_df