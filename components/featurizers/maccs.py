from src.chemical_representations import get_maccs_keys_fingerprint
import pandas as pd
import argparse
import numpy as np
from sklearn.feature_selection import VarianceThreshold


def main(args):
    """
    Compute MACCS fingerprints for SMILES strings in training and testing CSV files,
    remove low-variance features based on the training set, and save the results to new CSV files.
    """

    train_df = pd.read_csv(args.train_csv)
    train_fingerprints = []

    for smiles in train_df[args.smiles_col]:
        train_fp = get_maccs_keys_fingerprint(smiles)
        train_fingerprints.append(np.array(train_fp))

    fp_df = pd.DataFrame(
        train_fingerprints, columns=[f"maccs_key_{i}" for i in range(167)]
    )

    train_df = pd.concat([train_df, fp_df], axis=1)
    train_df = train_df.drop(columns=[args.smiles_col])

    test_df = pd.read_csv(args.test_csv)
    test_fingerprints = []

    for smiles in test_df[args.smiles_col]:
        test_fp = get_maccs_keys_fingerprint(smiles)
        test_fingerprints.append(np.array(test_fp))

    test_fp_df = pd.DataFrame(
        test_fingerprints, columns=[f"maccs_key_{i}" for i in range(167)]
    )

    test_df = pd.concat([test_df, test_fp_df], axis=1)
    test_df = test_df.drop(columns=[args.smiles_col])

    selector = VarianceThreshold(threshold=0.0)
    variance_check_df = train_df.drop(columns=args.activity_col)
    selector.fit(variance_check_df)
    selected_columns = variance_check_df.columns[selector.get_support()]

    train_df = train_df.drop(
        columns=variance_check_df.columns.difference(selected_columns)
    )
    test_df = test_df.drop(columns=test_df.columns.difference(selected_columns))

    train_df[args.activity_col] = pd.read_csv(args.train_csv)[args.activity_col]
    test_df[args.activity_col] = pd.read_csv(args.test_csv)[args.activity_col]

    train_df.to_csv(args.train_output_csv, index=False)
    test_df.to_csv(args.test_output_csv, index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="MACCS Keys Fingerprint Featurizer")

    parser.add_argument("--train_csv", type=str, help="Path to the training CSV file.")
    parser.add_argument("--test_csv", type=str, help="Path to the testing CSV file.")
    parser.add_argument(
        "--train_output_csv",
        type=str,
        help="Path to save the output training CSV file.",
    )
    parser.add_argument(
        "--test_output_csv", type=str, help="Path to save the output testing CSV file."
    )
    parser.add_argument(
        "--smiles_col",
        type=str,
        default="standardized_smiles",
        help="Column name for SMILES strings.",
    )
    parser.add_argument(
        "--activity_col",
        type=str,
        default=["antimicrobial_activity"],
        help="Column names for target variables.",
    )
    args = parser.parse_args()
    main(args)