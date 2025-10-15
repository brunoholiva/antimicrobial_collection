import pandas as pd
from src.utils.smiles_utils import SMILESValidator, Standardizer
from loguru import logger


class DataFrameSmilesProcessor:
    """
    Class for processing a DataFrame containing SMILES strings.
    """
    def __init__(self):
        self.validator = SMILESValidator()
        self.standardizer = Standardizer()

    def process(self, df:pd.DataFrame, smiles_column:str) -> pd.DataFrame:
        processed_df = df.copy()
        
        # 1. Validate initial SMILES
        logger.info(f"Starting SMILES validation of {len(processed_df)} entries.")

        processed_df['is_valid'] = processed_df[smiles_column].apply(self.validator.is_valid)
        num_invalid = len(processed_df) - processed_df['is_valid'].sum()

        logger.info(f"Number of invalid SMILES: {num_invalid}")
       
        valid_df = processed_df[processed_df['is_valid']].copy()
        
        logger.success(f"SMILES validation completed. {len(valid_df)} valid entries retained.")
    
        # 2. Standardize SMILES
        logger.info(f"Starting SMILES standardization of {len(valid_df)} entries.")

        valid_df['standardized_smiles'] = valid_df[smiles_column].apply(self.standardizer.standardize_smiles)
        valid_df['is_valid'] = valid_df['standardized_smiles'].apply(self.validator.is_valid)

        num_invalid_after_std = len(valid_df) - valid_df['is_valid'].sum()

        logger.info(f"Number of invalid SMILES after standardization: {num_invalid_after_std}")

        valid_df = valid_df[valid_df['is_valid']].copy()

        logger.success(f"SMILES standardization completed. {len(valid_df)} valid entries retained.")

        # 3. Remove duplicates based on standardized SMILES
        # TODO: Consider keeping the most relevant entry based on some criteria

        logger.info("Removing duplicates based on standardized SMILES.")
        valid_df = valid_df.drop_duplicates(subset='standardized_smiles', keep='first')
        logger.success(f"Duplicate removal completed. {len(valid_df)} unique entries retained, from {len(processed_df)} original entries.")

        cols_to_keep = ['antimicrobial_activity', 'standardized_smiles', 'target','source']
        valid_df = valid_df[[col for col in cols_to_keep if col in valid_df.columns]]

        return valid_df
    
def save_processed_df(df, name):
    df.to_csv(f"data/processed/{name}_processed.csv", index=False)
