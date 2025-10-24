import pandas as pd
from loguru import logger


def deduplicate_smiles(
    df: pd.DataFrame, rule, smiles_col: str = "standardized_smiles"
):
    """
    Deduplicate a DataFrame based on a SMILES column using a custom rule.

    Args:
        df: Input DataFrame containing SMILES strings.
        rule: Function defining the deduplication rule.
        smiles_col: Name of the column containing SMILES strings.

    Returns:
        A deduplicated DataFrame.
    """

    deduped_df = df.groupby(smiles_col, group_keys=False).apply(rule)
    return deduped_df.reset_index(drop=True)


def rule(group):
    """ 
    Custom deduplication rule for selecting a representative 
    row from a group of duplicates, when merging datasets.

    Args:
        group: DataFrame group containing duplicate SMILES entries.

    Returns:
        A single row (as a Series) representing the group.
    """
    # 1. Prefer co_add
    co_add_rows = group[group["source"] == "CO_ADD"]
    if not co_add_rows.empty:
        return co_add_rows.iloc[0]

    # 2. Prefer datasets with absolute cutoff
    absolute_cutoff_sources = {"stokes_2020", "wong_2024"}
    abs_cutoff_rows = group[group["source"].isin(absolute_cutoff_sources)]
    if not abs_cutoff_rows.empty:
        return abs_cutoff_rows.iloc[0]

    # 3. Prefer the largest dataset (by source frequency in the whole dataset)
    source_counts = group["source"].value_counts()
    largest_source = source_counts.idxmax()
    largest_rows = group[group["source"] == largest_source]
    return largest_rows.iloc[0]


def merge_and_deduplicate(dfs, smiles_col, rule):
    merged = pd.concat(dfs, ignore_index=True)
    deduped = deduplicate_smiles(merged, rule, smiles_col)
    logger.success(
        f"Merged dataset has {len(merged)} entries; after deduplication, {len(deduped)} unique entries remain. (Removed {len(merged) - len(deduped)} duplicates.)"
    )
    return deduped
