from rdkit import Chem
import pandas as pd
from rdkit.Chem.MolStandardize import rdMolStandardize

class SMILESValidator:
    """
    Class for validating SMILES strings.
    """
    def is_valid(self, smiles:str) -> bool:
        """
        Check if a SMILES string is valid by attempting to create an RDKit molecule object from it
        Args:
            smiles (str): The SMILES string to validate.
        Returns:
            True if the SMILES string is valid, False otherwise.
        """

        if not isinstance(smiles, str):
            return False
        molecule = Chem.MolFromSmiles(smiles)
        return molecule is not None
    

class Standardizer:
    """
    Class for standardizing molecules.
    """
    def __init__(self):
        self.uncharger = rdMolStandardize.Uncharger()
        self.tautomer_enumerator = rdMolStandardize.TautomerEnumerator()


    def standardize_mol(self, mol:Chem.Mol) -> Chem.Mol:
        """
        Standardize a molecule using RDKit's MolStandardize module.
        """
        if mol is None:
            return None

        clean_molecule = rdMolStandardize.Cleanup(mol)
        parent_clean_molecule = rdMolStandardize.FragmentParent(clean_molecule)
        uncharged_clean_molecule = self.uncharger.uncharge(parent_clean_molecule)
        tautomer_uncharged_clean_molecule = self.tautomer_enumerator.Canonicalize(uncharged_clean_molecule)

        return tautomer_uncharged_clean_molecule
    

    def standardize_smiles(self, smiles:str) -> str:
        """
        Standardize a SMILES string.
        """
        molecule = Chem.MolFromSmiles(smiles)
        if molecule:
            standardized_molecule = self.standardize_mol(molecule)
            return Chem.MolToSmiles(standardized_molecule)
        return None