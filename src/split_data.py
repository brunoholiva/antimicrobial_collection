# this script provides ways to split chemical data into training and test sets.
from sklearn.model_selection import train_test_split
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit import DataStructs
from rdkit.ML.Cluster import Butina
from src.chemical_representations import get_morgan_fingerprint
import deepchem as dc
from scipy.spatial.distance import pdist


def random_split(
    df: pd.DataFrame,
    test_size: float = 0.2,
    random_state: int = 333,
    stratify_col: str = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Splits a DataFrame into training and test sets using random sampling.

    Args:
        df: A dataframe containing the data to be split.
        test_size: The proportion of the dataset to include in the test split.
        random_state: Random seed for reproducibility.
        stratify_col: Column name to use for stratification. If None, no stratification is applied.
    """

    return train_test_split(
        df,
        test_size=test_size,
        random_state=random_state,
        shuffle=True,
        stratify=df[stratify_col] if stratify_col else None,
    )


def generate_murcko_scaffold(smiles: str) -> str:
    """
    Generates a Bemis-Murcko scaffold from a SMILES string.

    Args:
        smiles: A SMILES string representing a molecule.

    Returns:
        A SMILES string representing the Bemis-Murcko scaffold of the molecule.
        If the scaffold cannot be generated, returns the original SMILES string.
    """

    mol = Chem.MolFromSmiles(smiles)
    if mol:
        scaffold = MurckoScaffold.MurckoScaffoldSmiles(mol=mol, includeChirality=False)
        if len(scaffold) == 0:
            scaffold = smiles
    return scaffold


def murcko_scaffold_split(
    df: pd.DataFrame,
    smiles_col: str = "standardized_smiles",
    activity_col: str = "antimicrobial_activity",
    test_size: float = 0.2,
    random_state: int = 333,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Splits a DataFrame into training and test sets based on Bemis-Murcko scaffolds.

    Args:
        df: A dataframe containing the data to be split.
        smiles_col: The name of the column containing SMILES strings.
        activity_col: The name of the column containing activity labels (0 or 1).
        test_size: The proportion of the dataset to include in the test split.
        random_state: Random seed for reproducibility.

    Returns:
        train_df: Training set DataFrame.
        test_df: Test set DataFrame.
    """

    df["scaffold"] = df[smiles_col].apply(generate_murcko_scaffold)

    scaffold_groups = (
        df.groupby("scaffold")
        .agg(has_active=(activity_col, "max"), mol_count=("scaffold", "size"))
        .reset_index()
    )

    train_scaffolds, test_scaffolds = train_test_split(
        scaffold_groups,
        test_size=test_size,
        random_state=random_state,
        shuffle=True,
        # Stratify by whether the scaffold contains an active
        stratify=scaffold_groups["has_active"],
    )

    train_scaffold_list = train_scaffolds["scaffold"]
    test_scaffold_list = test_scaffolds["scaffold"]

    train_df = df[df["scaffold"].isin(train_scaffold_list)].drop(columns=["scaffold"])
    test_df = df[df["scaffold"].isin(test_scaffold_list)].drop(columns=["scaffold"])

    return train_df, test_df


def calculate_tanimoto_distances(fingerprint_list: list) -> np.ndarray:
    """
    Calculates the Tanimoto distance matrix using scipy's pdist.

    Args:
        fingerprint_list: A list of RDKit ExplicitBitVect objects.

    Returns:
        A 1D numpy array containing the condensed pairwise Tanimoto distances.
    """
    if not fingerprint_list:
        return np.array([])

    np_fps = np.array([list(fp) for fp in fingerprint_list], dtype=bool)

    # Calculate pairwise Jaccard distances (equivalent to Tanimoto distance)
    return pdist(np_fps, metric="jaccard")


def butina_clustering(
    fingerprint_list: list[Chem.DataStructs.ExplicitBitVect], cutoff: float = 0.65
) -> list[int]:
    """
    Clusters Morgan fingerprints using the Butina algorithm.

    Args:
        fingerprint_list: List of RDKit fingerprint objects.
        cutoff: Tanimoto distance cutoff for clustering.

    Returns:
        cluster_id_list: List of cluster IDs corresponding to each fingerprint.
    """

    distances = calculate_tanimoto_distances(fingerprint_list)
    n_fingerprints = len(fingerprint_list)

    clusters = Butina.ClusterData(distances, n_fingerprints, cutoff, isDistData=True)

    cluster_id_list = np.zeros(n_fingerprints, dtype=int)
    for cluster_num, cluster in enumerate(clusters):
        for member in cluster:
            cluster_id_list[member] = cluster_num

    return cluster_id_list.tolist()


def butina_split(
    df: pd.DataFrame,
    smiles_col: str = "standardized_smiles",
    activity_col: str = "antimicrobial_activity",
    test_size: float = 0.2,
    random_state: int = 333,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Splits a DataFrame into training and test sets based on Butina clustering of Morgan fingerprints.

    Args:
        df: A dataframe containing the data to be split.
        smiles_col: The name of the column containing SMILES strings.
        activity_col: The name of the column containing activity labels (0 or 1).
        test_size: The proportion of the dataset to include in the test split.
        random_state: Random seed for reproducibility.

    Returns:
        train_df: Training set DataFrame.
        test_df: Test set DataFrame.

    """

    fingerprint_list = [get_morgan_fingerprint(s) for s in df[smiles_col]]
    df["cluster_id"] = butina_clustering(fingerprint_list, cutoff=0.65)

    cluster_groups = (
        df.groupby("cluster_id")
        .agg(has_active=(activity_col, "max"), mol_count=("cluster_id", "size"))
        .reset_index()
    )
    train_scaffolds, test_scaffolds = train_test_split(
        cluster_groups,
        test_size=test_size,
        random_state=random_state,
        shuffle=True,
        # Stratify by whether the scaffold contains an active
        stratify=cluster_groups["has_active"],
    )

    train_cluster_list = train_scaffolds["cluster_id"]
    test_cluster_list = test_scaffolds["cluster_id"]

    train_df = df[df["cluster_id"].isin(train_cluster_list)].drop(
        columns=["cluster_id"]
    )
    test_df = df[df["cluster_id"].isin(test_cluster_list)].drop(columns=["cluster_id"])

    return train_df, test_df


def maxmin_split(
    df: pd.DataFrame,
    smiles_col: str = "standardized_smiles",
    activity_col: str = "antimicrobial_activity",
    train_size: float = 0.8,
    random_state: int = 333,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Splits a DataFrame into training and test sets using DeepChem"s MaxMinSplitter.

    Args:
        df: A dataframe containing the data to be split.
        smiles_col: The name of the column containing SMILES strings.
        activity_col: The name of the column containing activity labels (0 or 1).
        train_size: The proportion of the dataset to include in the train split.
        random_state: Random seed for reproducibility.

    Returns:
        train_df: Training set DataFrame.
        test_df: Test set DataFrame.
    """

    smiles = df[smiles_col].astype(str).values
    labels = df[activity_col].values
    X_dummy = np.zeros((len(df), 1))

    dataset = dc.data.NumpyDataset(X=X_dummy, y=labels, ids=smiles)

    maxminsplitter = dc.splits.MaxMinSplitter()
    train_dataset, test_dataset = maxminsplitter.train_test_split(
        dataset, frac_train=train_size, random_state=random_state
    )

    train_df = (
        df[df["standardized_smiles"].astype(str).isin(train_dataset.ids)]
        .reset_index(drop=True)
        .drop(columns=["scaffold", "cluster_id"])
    )
    test_df = (
        df[df["standardized_smiles"].astype(str).isin(test_dataset.ids)]
        .reset_index(drop=True)
        .drop(columns=["scaffold", "cluster_id"])
    )

    return train_df, test_df
