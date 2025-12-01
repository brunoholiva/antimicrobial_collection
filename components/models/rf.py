import argparse
from sklearn.ensemble import RandomForestClassifier
from sklearn.pipeline import Pipeline
from sklearn.model_selection import PredefinedSplit
from sklearn.feature_selection import SelectKBest, VarianceThreshold, f_classif
from skopt import BayesSearchCV
from skopt.space import Integer, Categorical
import pandas as pd
import joblib


def main(args):
    df = pd.read_csv(args.input_csv)
    
    split_index = df["cv_fold"].values
    cv_strategy = PredefinedSplit(test_fold=split_index)
    X = df.drop(columns=[args.activity_col, "cv_fold"])
    y = df[args.activity_col]

    pipeline = Pipeline(
        [
            ("variance", VarianceThreshold()),
            ("select", SelectKBest(score_func=f_classif)),
            ("model",RandomForestClassifier(
                    random_state=args.random_state, 
                    class_weight="balanced", 
                    n_jobs=args.n_jobs
                ),
            ),
        ]
    )

    X_var = VarianceThreshold().fit_transform(X)

    opt = BayesSearchCV(
        estimator=pipeline,
        search_spaces={
            "select__k": Integer(16, X_var.shape[1]),
            "model__n_estimators": Integer(100, 1000),
            "model__max_depth": Integer(5, 50),
            "model__min_samples_split": Integer(2, 10),
            "model__min_samples_leaf": Integer(1, 4),
            "model__bootstrap": Categorical([True, False]),
            "model__max_features": Categorical(["sqrt", "log2"]),
        },
        scoring="average_precision",
        cv=cv_strategy,
        n_iter=args.n_iter,
        n_jobs=1,
        random_state=args.random_state,
        verbose=2
    )
    opt.fit(X, y)
    model = opt.best_estimator_
    joblib.dump(model, args.output_model_path)


    cv_df = pd.DataFrame(opt.cv_results_)
    cols_to_keep = [col for col in cv_df.columns if "mean_test" in col or "std_test" in col or "params" in col]
    cv_df_clean = cv_df[cols_to_keep].copy()
    cv_df_clean['is_best'] = False
    cv_df_clean.loc[opt.best_index_, 'is_best'] = True
    cv_df_clean.to_csv(args.output_metrics, index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Random Forest Classifier Trainer")
    parser.add_argument("--input_csv", type=str, help="Path to the training CSV file.")
    parser.add_argument("--output_model_path", type=str, help="Path to save the trained model.")
    parser.add_argument("--activity_col", type=str, default="antimicrobial_activity", help="Column name for activity labels.")
    parser.add_argument("--output_metrics", type=str, required=True)
    parser.add_argument("--random_state", type=int, default=333, help="Random seed for reproducibility.")
    parser.add_argument("--n_jobs", type=int, help="Number of parallel jobs to run.")
    parser.add_argument("--n_iter", type=int, default=500, help="Number of iterations for Bayesian optimization.")
    args = parser.parse_args()
    main(args)