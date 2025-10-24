from src.split_data import butina_split
import argparse
import pandas as pd


def main(args):
    data = pd.read_csv(args.input_csv)
    train_set, test_set = butina_split(
        data,
        activity_col=args.activity_col,
        smiles_col=args.smiles_col,
        random_state=args.random_state,
        test_size=args.test_size,
    )
    train_set = train_set[[args.smiles_col, args.activity_col]]
    test_set = test_set[[args.smiles_col, args.activity_col]]
    
    train_set.to_csv(args.train_output_csv, index=False)
    test_set.to_csv(args.test_output_csv, index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Split dataset using Taylor-Butina algorithm"
    )
    parser.add_argument(
        "--input_csv", type=str, required=True, help="Path to input CSV file"
    )
    parser.add_argument(
        "--train_output_csv",
        type=str,
        required=True,
        help="Path to train output CSV file",
    )
    parser.add_argument(
        "--test_output_csv",
        type=str,
        required=True,
        help="Path to test output CSV file",
    )
    parser.add_argument(
        "--activity_col",
        type=str,
        default="antimicrobial_activity",
        help="Column name for activity labels",
    )
    parser.add_argument(
        "--smiles_col",
        type=str,
        default="standardized_smiles",
        help="Column name for SMILES strings",
    )
    parser.add_argument(
        "--test_size",
        type=float,
        default=0.2,
        help="Proportion of the dataset to include in the test split",
    )
    parser.add_argument(
        "--random_state", type=int, default=333, help="Random seed for reproducibility"
    )
    args = parser.parse_args()
    main(args)
