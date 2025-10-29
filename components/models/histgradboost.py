from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.model_selection import RandomizedSearchCV
import pandas as pd
import joblib
import argparse


def main(args):
    train_df = pd.read_csv(args.train_csv)

    X_train = train_df.drop(columns=[args.activity_col])
    y_train = train_df[args.activity_col]

    model = HistGradientBoostingClassifier(
        random_state=args.random_state, class_weight="balanced"
    )

    param_grid = {
        "learning_rate": [0.01, 0.1, 0.2, 0.3],
        "max_iter": [100, 200, 300, 400, 500],
        "min_samples_leaf": [20, 50, 100, 200],
        "l2_regularization": [0.0, 0.1, 0.2, 0.3],
        "min_samples_leaf": [15, 20, 25, 30, 40],
    }

    grid_search = RandomizedSearchCV(
        estimator=model,
        param_distributions=param_grid,
        cv=10,
        scoring="average_precision",
        n_iter=500,
        n_jobs=args.n_jobs,
        verbose=2,
    )

    grid_search.fit(X_train, y_train)
    model = grid_search.best_estimator_
    model.fit(X_train, y_train)
    joblib.dump(model, args.output_model_path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Histogram-based Gradient Boosting Classifier Trainer"
    )
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
        default=333,
        help="Random state for reproducibility.",
    )
    parser.add_argument(
        "--n_jobs",
        type=int,
        default=1,
        help="Number of parallel jobs to run.",
    )
    args = parser.parse_args()
    main(args)
