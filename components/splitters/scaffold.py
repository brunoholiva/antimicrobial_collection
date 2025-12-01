import argparse
import pandas as pd
from rdkit.Chem.Scaffolds import MurckoScaffold
from sklearn.model_selection import StratifiedGroupKFold


def main(args):
    df = pd.read_csv(args.input_csv)

    scaffolds = df[args.smiles_col].apply(
        lambda x: MurckoScaffold.MurckoScaffoldSmiles(mol=None, smiles=x)
    )

    gkf = StratifiedGroupKFold(n_splits=args.k_folds, shuffle=True, random_state=args.random_state)
    df["cv_fold"] = -1

    X_dummy = df.drop(columns=[args.activity_col])
    y_dummy = df[args.activity_col]

    for fold_id, (train_idx, val_idx) in enumerate(
        gkf.split(X_dummy, y_dummy, groups=scaffolds)
    ):
        df.iloc[val_idx, df.columns.get_loc("cv_fold")] = fold_id

    df.to_csv(args.output_csv, index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Scaffold-based K-Fold Splitter")
    parser.add_argument("--input_csv", type=str, required=True, help="Path to input CSV file")
    parser.add_argument("--output_csv", type=str, required=True, help="Path to output CSV file")
    parser.add_argument("--smiles_col", type=str, default="standardized_smiles", help="Column name for SMILES strings")
    parser.add_argument("--activity_col", type=str, default="antimicrobial_activity", help="Column name for activity labels")
    parser.add_argument("--random_state", type=int, default=333, help="Random seed for reproducibility")
    parser.add_argument("--k_folds", type=int, default=3, help="Number of folds for K-Fold splitting")
    args = parser.parse_args()
    main(args)