from src.chemical_representations import get_morgan_fingerprint
import pandas as pd
import argparse
import numpy as np


def main(args):
    """
    Compute Morgan fingerprints for SMILES strings in a CSV file (dataframe) and save the results to a new CSV file.

    Args:
        args: Command line arguments containing input and output file paths, SMILES column name,
              fingerprint radius, and number of bits.
    """
    df = pd.read_csv(args.input_csv)
    fingerprints = []
    for smiles in df[args.smiles_col]:
        fp = get_morgan_fingerprint(
            smiles,
            radius=args.radius,
            n_bits=args.n_bits,
        )
        fingerprints.append(np.array(fp))

    fp_df = pd.DataFrame(
        fingerprints, columns=[f"morgan_bit_{i}" for i in range(args.n_bits)]
    )
    df = pd.concat([df, fp_df], axis=1)
    df = df.drop(columns=[args.smiles_col])
    df.to_csv(args.output_csv, index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Morgan Fingerprint Featurizer")
    parser.add_argument("--input_csv", type=str, help="Path to the input CSV file.")
    parser.add_argument(
        "--output_csv", type=str, help="Path to save the output CSV file."
    )
    parser.add_argument(
        "--smiles_col",
        type=str,
        default="standardized_smiles",
        help="Column name for SMILES strings.",
    )
    parser.add_argument(
        "--radius",
        type=int,
        default=2,
        help="Radius for Morgan fingerprint.",
    )
    parser.add_argument(
        "--n_bits",
        type=int,
        default=2048,
        help="Number of bits for Morgan fingerprint.",
    )

    args = parser.parse_args()
    main(args)
