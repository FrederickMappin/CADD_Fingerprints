from abc import ABC, abstractmethod
import pandas as pd

__all__ = ['Similarity']

class Similarity(ABC):
    @abstractmethod
    def calculate_similarity(self, fingerprints_df):
        """
        Calculate the similarity for the given fingerprints DataFrame.
        
        Parameters:
        fingerprints_df (pd.DataFrame): DataFrame containing the fingerprints.
        
        Returns:
        pd.DataFrame: DataFrame containing the similarity scores.
        """
        pass