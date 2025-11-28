from sklearn.gaussian_process import GaussianProcessClassifier
from sklearn.feature_selection import SelectKBest, VarianceThreshold, f_classif
from skopt import BayesSearchCV
from skopt.space import Real, Integer, Categorical
from sklearn.pipeline import Pipeline
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
            ("model", GaussianProcessClassifier(random_state=args.random_state, n_jobs=args.n_jobs)),
        ]
    )

    opt = BayesSearchCV(
        estimator=pipeline,
        search_spaces={
            "select__k": Integer(16, X_train.shape[1]),
            "model__optimizer": Categorical(["fmin_l_bfgs_b", None]),
            "model__max_iter_predict": Integer(50, 400),
        },
        cv=3,
        n_iter=50,
        scoring="average_precision",
        n_jobs=1,
        verbose=2,
        random_state=args.random_state
    )


    opt.fit(X_train, y_train)
    model = opt.best_estimator_
    joblib.dump(model, args.output_model_path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Gaussian Process Classifier Trainer")
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
