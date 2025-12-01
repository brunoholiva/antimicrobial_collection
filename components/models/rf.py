import argparse
from sklearn.ensemble import RandomForestClassifier
from sklearn.pipeline import Pipeline
from sklearn.feature_selection import SelectKBest, VarianceThreshold, f_classif
from skopt import BayesSearchCV
from skopt.space import Integer, Categorical
import pandas as pd
import joblib
import numpy as np


def main(args):
    train_df = pd.read_csv(args.train_csv)

    X_train = train_df.drop(columns=[args.activity_col])
    y_train = train_df[args.activity_col]

    pipeline = Pipeline(
        [
            ("variance", VarianceThreshold()),
            ("select", SelectKBest(score_func=f_classif)),
            (
                "model",
                RandomForestClassifier(
                    random_state=args.random_state, class_weight="balanced", n_jobs=args.n_jobs
                ),
            ),
        ]
    )

    opt = BayesSearchCV(
        estimator=pipeline,
        search_spaces={
            "select__k": Integer(16, X_train.shape[1]),
            "model__n_estimators": Integer(100, 1000),
            "model__max_depth": Integer(5, 50),
            "model__min_samples_split": Integer(2, 10),
            "model__min_samples_leaf": Integer(1, 4),
            "model__bootstrap": Categorical([True, False]),
            "model__max_features": Categorical(["sqrt", "log2"]),
        },
        scoring="average_precision",
        cv=3,
        n_iter=300,
        n_jobs=1,
        random_state=args.random_state
    )

    opt.fit(X_train, y_train)
    model = opt.best_estimator_
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
