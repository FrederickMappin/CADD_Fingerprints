from abc import ABC, abstractmethod
import pandas as pd

__all__ = ['Similarity']

class Similarity(ABC):
    @abstractmethod
    def calculate_similarity(self, input_file, output_file):
        """
        Calculate the similarity for the given fingerprints CSV file and save the result to the output file.
        
        Parameters:
        input_file (str): Path to the input file containing fingerprints.
        output_file (str): Path to the output file where the similarity matrix will be saved.
        
        Returns:
        pd.DataFrame: DataFrame containing the similarity scores.
        """
        pass