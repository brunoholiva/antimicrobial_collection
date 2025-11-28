from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.feature_selection import SelectKBest, VarianceThreshold, f_classif
from sklearn.pipeline import Pipeline
from skopt import BayesSearchCV
from skopt.space import Real, Integer
import pandas as pd
import joblib
import argparse


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
                HistGradientBoostingClassifier(
                    random_state=args.random_state, class_weight="balanced"
                ),
            ),
        ]
    )

    opt = BayesSearchCV(
        estimator=pipeline,
        search_spaces={
            "select__k": Integer(16, X_train.shape[1]),
            "model__learning_rate": Real(0.01, 0.3, prior="log-uniform"),
            "model__max_iter": Integer(100, 500),
            "model__min_samples_leaf": Integer(20, 200),
            "model__l2_regularization": Real(0.0, 0.3),
        },
        cv=3,
        n_iter=300,
        scoring="average_precision",
        n_jobs=args.n_jobs,
        verbose=2,
        random_state=args.random_state
    )

    opt.fit(X_train, y_train)
    model = opt.best_estimator_
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
