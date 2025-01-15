from abc import ABC, abstractmethod

__all__ = ['Fingerprint']

class Fingerprint(ABC):
    @abstractmethod
    def calculate(self, input_file, output_file):
        """
        Calculate the fingerprint for the given input file and save the result to the output file.
        
        Parameters:
        input_file (str): Path to the input file containing SMILES strings.
        output_file (str): Path to the output file where the fingerprints will be saved.
        
        Returns:
        pd.DataFrame: DataFrame containing the original data and the calculated fingerprints.
        """
        pass