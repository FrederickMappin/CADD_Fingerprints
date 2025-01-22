import pandas as pd
from rdkit import DataStructs
from rdkit.DataStructs import ExplicitBitVect
from similarity.similarity import Similarity

class TanimotoSimilarity(Similarity):
    def calculate_similarity(self, input_file, output_file):
        # Read the input CSV file containing fingerprints
        fingerprints_df = pd.read_csv(input_file)
        
        # Check if 'molecule_chembl_id' column exists
        if 'molecule_chembl_id' not in fingerprints_df.columns:
            raise ValueError("Input file must contain a 'molecule_chembl_id' column")
        
        # Set 'molecule_chembl_id' as the index
        fingerprints_df.set_index('molecule_chembl_id', inplace=True)
        
        # Extract fingerprints
        fingerprints = []
        for i, row in fingerprints_df.iterrows():
            fp = ExplicitBitVect(len(row))
            for idx, bit in enumerate(row):
                if bit == 1:
                    fp.SetBit(idx)
            fingerprints.append(fp)
        
        # Calculate Tanimoto similarity
        num_fingerprints = len(fingerprints)
        similarity_matrix = pd.DataFrame(index=fingerprints_df.index, columns=fingerprints_df.index)
        
        for i in range(num_fingerprints):
            for j in range(i, num_fingerprints):
                similarity = DataStructs.TanimotoSimilarity(fingerprints[i], fingerprints[j])
                similarity_matrix.iloc[i, j] = round(similarity, 2)
                similarity_matrix.iloc[j, i] = round(similarity, 2)
        
        # Save the similarity matrix to the output CSV file with two decimal places
        similarity_matrix.to_csv(output_file, float_format='%.2f')
        
        # Return the similarity matrix DataFrame
        return similarity_matrix

if __name__ == "__main__":
    import sys
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    ts = TanimotoSimilarity()
    similarity_matrix = ts.calculate_similarity(input_file, output_file)
    print(f"Similarity matrix saved to {output_file}")