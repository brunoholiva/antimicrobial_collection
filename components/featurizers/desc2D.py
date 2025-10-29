import pandas as pd
import argparse
import numpy as np
from sklearn.preprocessing import StandardScaler

# Assuming get_desc2d is in src.chemical_representations
from src.chemical_representations import get_desc2d


def main(args):
    """
    Compute and normalize 2D descriptors for SMILES strings from train and test CSVs.
    """
    train_df = pd.read_csv(args.train_csv)
    test_df = pd.read_csv(args.test_csv)

    # Remove activity column before featurization
    train_activity = train_df[args.activity_col]
    test_activity = test_df[args.activity_col]
    train_df = train_df.drop(columns=[args.activity_col])
    test_df = test_df.drop(columns=[args.activity_col])

    train_descriptors_list = []
    for smiles in train_df[args.smiles_col]:
        descriptor = get_desc2d(smiles)
        train_descriptors_list.append(descriptor)

    desc_cols = [f"desc2d_{i}" for i in range(train_descriptors_list[0].shape[0])]
    train_desc_df = pd.DataFrame(train_descriptors_list, columns=desc_cols)

    test_descriptors_list = []
    for smiles in test_df[args.smiles_col]:
        descriptor = get_desc2d(smiles)
        test_descriptors_list.append(descriptor)
    test_desc_df = pd.DataFrame(test_descriptors_list, columns=desc_cols)

    # drop any columns with NaNs
    train_desc_df = train_desc_df.dropna(axis=1)
    test_desc_df = test_desc_df[train_desc_df.columns]
    test_desc_df = test_desc_df.dropna(axis=1)
    train_desc_df = train_desc_df[test_desc_df.columns]

    desc_cols = train_desc_df.columns

    scaler = StandardScaler()
    train_desc_scaled = scaler.fit_transform(train_desc_df)
    test_desc_scaled = scaler.transform(test_desc_df)

    train_desc_scaled_df = pd.DataFrame(train_desc_scaled, columns=desc_cols)
    test_desc_scaled_df = pd.DataFrame(test_desc_scaled, columns=desc_cols)

    train_output_df = pd.concat([train_desc_scaled_df, train_activity.reset_index(drop=True)], axis=1)
    train_output_df.to_csv(args.train_output_csv, index=False)

    test_output_df = pd.concat([test_desc_scaled_df, test_activity.reset_index(drop=True)], axis=1)
    test_output_df.to_csv(args.test_output_csv, index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="2D Descriptor Featurizer with Normalization"
    )

    parser.add_argument(
        "--train_csv", type=str, help="Path to the input training CSV file."
    )
    parser.add_argument("--test_csv", type=str, help="Path to the input test CSV file.")
    parser.add_argument(
        "--train_output_csv",
        type=str,
        help="Path to save the output training features.",
    )
    parser.add_argument(
        "--test_output_csv", type=str, help="Path to save the output test features."
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
        default="antimicrobial_activity",
        help="Column name for activity values.",
    )
    args = parser.parse_args()
    main(args)
