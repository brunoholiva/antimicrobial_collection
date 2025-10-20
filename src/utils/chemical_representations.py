from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator, MACCSkeys


def get_morgan_fingerprint(
    smiles: str, radius: int = 2, n_bits: int = 2048
) -> Chem.DataStructs.ExplicitBitVect:
    """
    Generate a Morgan fingerprint from a SMILES string.

    Args:
        smiles: SMILES string.
        radius: The radius around each atom for the Morgan Generator.
        n_bits: Total number of bits of the fingerprint.

    Returns:
        fingerprint: A Fingerprint (ExplicitBitVect)
    """

    morgan_generator = rdFingerprintGenerator.GetMorganGenerator(radius=2)
    molecule = Chem.MolFromSmiles(smiles)
    if molecule:
        fingerprint = morgan_generator.GetFingerprint(molecule)
    return fingerprint


def get_maccs_keys_fingerprint(smiles: str) -> Chem.DataStructs.ExplicitBitVect:
    """
    Generate a MACCS keys fingerprint from a SMILES string.
    Args:
        smiles: SMILES string.
    Returns:
        fingerprint: A Fingerprint (ExplicitBitVect)
    """

    molecule = Chem.MolFromSmiles(smiles)
    if molecule:
        fingerprint = MACCSkeys.GenMACCSKeys(molecule)
    return fingerprint
