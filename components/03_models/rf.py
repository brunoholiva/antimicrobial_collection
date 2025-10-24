import argparse
from sklearn.ensemble import RandomForestClassifier
import pandas as pd
import joblib


def main(args):
    train_df = pd.read_csv(args.train_csv)

    X_train = train_df.drop(columns=[args.activity_col])
    y_train = train_df[args.activity_col]

    model = RandomForestClassifier(
        random_state=args.random_state,
        n_jobs=args.n_jobs,
        n_estimators=args.n_estimators,
        class_weight="balanced"
    )

    model.fit(X_train, y_train)

    joblib.dump(model, args.output_model_path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Random Forest Classifier Trainer")
    parser.add_argument("--train_csv", type=str, help="Path to the training CSV file.")
    parser.add_argument(
        "--output_model_path", type=str, help="Path to save the trained model."
    )
    parser.add_argument(
        "--activity_col",
        type=str,
        default="antimicrobial_activity",
        help="Column name for activity labels.",
    )
    parser.add_argument(
        "--random_state",
        type=int,
        default=42,
        help="Random seed for reproducibility.",
    )
    parser.add_argument(
        "--n_jobs",
        type=int,
        default=-1,
        help="Number of parallel jobs to run.",
    )
    parser.add_argument(
        "--n_estimators",
        type=int,
        default=100,
        help="Number of trees in the random forest.",
    )

    args = parser.parse_args()
    main(args)
