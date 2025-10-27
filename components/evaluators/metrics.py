import argparse
import pandas as pd
import joblib
from sklearn.metrics import roc_auc_score, accuracy_score, precision_score, recall_score, average_precision_score


def main(args):
    model = joblib.load(args.model_path)
    test_df = pd.read_csv(args.test_csv)

    X_test = test_df.drop(columns=[args.activity_col])
    y_test = test_df[args.activity_col]

    y_pred_proba = model.predict_proba(X_test)[:, 1]
    y_pred_class = model.predict(X_test)

    metrics = {
        "test_auc": roc_auc_score(y_test, y_pred_proba),
        "test_accuracy": accuracy_score(y_test, y_pred_class),
        "test_precision": precision_score(y_test, y_pred_class),
        "test_recall": recall_score(y_test, y_pred_class),
        "test_average_precision": average_precision_score(y_test, y_pred_proba)
    }

    results_data = {
        "dataset": args.dataset_name,
        "split": args.split_name,
        "featurizer": args.featurizer_name,
        "model": args.model_name,
        "evaluator": args.evaluator_name,
        "random_state": args.random_state,
        **metrics,
    }

    results_df = pd.DataFrame([results_data])
    results_df.to_csv(args.output_results, index=False)
    print(f"Results saved to {args.output_results}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Model Evaluation Metrics Calculator")
    parser.add_argument(
        "--model_path", type=str, help="Path to the trained model file."
    )
    parser.add_argument("--test_csv", type=str, help="Path to the test CSV file.")
    parser.add_argument(
        "--activity_col",
        type=str,
        default="antimicrobial_activity",
        help="Column name for activity labels.",
    )
    parser.add_argument("--dataset_name", type=str, help="Name of the dataset.")
    parser.add_argument("--split_name", type=str, help="Name of the data split method.")
    parser.add_argument(
        "--featurizer_name", type=str, help="Name of the featurizer used."
    )
    parser.add_argument("--model_name", type=str, help="Name of the model used.")
    parser.add_argument(
        "--evaluator_name", type=str, help="Name of the evaluator used."
    )
    parser.add_argument("--random_state", type=int, default=42, help="Random seed used.")
    parser.add_argument(
        "--output_results", type=str, help="Path to save the evaluation results CSV."
    )

    args = parser.parse_args()
    main(args)
