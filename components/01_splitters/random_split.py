from src.split_data import random_split
import argparse
import pandas as pd


def main(args):
    df = pd.read_csv(args.input_csv)
    train_df, test_df = random_split(
        df,
        test_size=args.test_size,
        random_state=args.random_state,
        activity_col=args.activity_col,
        smiles_col=args.smiles_col
    )
    train_df = train_df[[args.smiles_col, args.activity_col]]
    test_df = test_df[[args.smiles_col, args.activity_col]]

    train_df.to_csv(args.train_output_csv, index=False)
    test_df.to_csv(args.test_output_csv, index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Murcko Scaffold Splitter")
    parser.add_argument("--input_csv", type=str, help="Path to the input CSV file.")
    parser.add_argument(
        "--train_output_csv", type=str, help="Path to save the training set CSV file."
    )
    parser.add_argument(
        "--test_output_csv", type=str, help="Path to save the test set CSV file."
    )
    parser.add_argument(
        "--activity_col",
        type=str,
        default="antimicrobial_activity",
        help="Column name for activity labels.",
    )
    parser.add_argument(
        "--smiles_col",
        type=str,
        default="standardized_smiles",
        help="Column name for SMILES strings.",
    )
    parser.add_argument(
        "--test_size",
        type=float,
        default=0.2,
        help="Proportion of the dataset to include in the test split.",
    )
    parser.add_argument(
        "--random_state", type=int, default=333, help="Random seed for reproducibility."
    )

    args = parser.parse_args()
    main(args)
