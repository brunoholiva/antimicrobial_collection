import argparse
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import RandomizedSearchCV
import pandas as pd
import joblib
import numpy as np

def main(args):
    train_df = pd.read_csv(args.train_csv)

    X_train = train_df.drop(columns=[args.activity_col])
    y_train = train_df[args.activity_col]

    model = RandomForestClassifier(
        random_state=args.random_state,
        class_weight="balanced_subsample"
    )
    
    param_grid = {
        "n_estimators": [int(x) for x in np.linspace(start = 200, stop = 2000, num = 10)],
        "max_features": ["sqrt", "log2"],
        "min_samples_split": [2, 5, 10],
        "min_samples_leaf": [1, 2, 4],
        "bootstrap": [True, False]
    }

    grid_search = RandomizedSearchCV(
        estimator=model,
        param_distributions=param_grid,
        cv=5,
        scoring="average_precision",
        n_iter=500,
        n_jobs=args.n_jobs,
        verbose=2
    )
    
    grid_search.fit(X_train, y_train)
    model = grid_search.best_estimator_
    
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
        default=333,
        help="Random seed for reproducibility.",
    )
    parser.add_argument(
        "--n_jobs",
        type=int,
        default=1,
        help="Number of parallel jobs to run.",
    )


    args = parser.parse_args()
    main(args)
