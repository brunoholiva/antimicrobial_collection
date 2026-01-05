import pandas as pd
import argparse
import ast

def main(args):
    df = pd.read_csv(args.input_csv)

    if 'is_best' in df.columns:
        best_row = df[df['is_best'] == True].iloc[0]
    else:
        best_row = df.loc[df['mean_test_score'].idxmax()]

    params_str = best_row['params'].replace("OrderedDict", "")
    params_dict = ast.literal_eval(params_str)    
    params_dict = dict(params_dict)

    
    clean_params = {}
    for k, v in params_dict.items():
        clean_key = k.replace("model__", "").replace("select__", "")
        clean_params[clean_key] = v

    final_data = {
        "Dataset": args.dataset_name,
        "Splitter": args.splitter_name,
        "Featurizer": args.featurizer_name,
        "Model": args.model_name,
        "Model_characteristics": args.model_path,
        "Mean_AP": round(best_row['mean_test_score'], 4),
        "Std_AP": round(best_row['std_test_score'], 4),
    }
    

    final_df = pd.DataFrame([final_data])
    final_df.to_csv(args.output_csv, index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_csv", type=str, required=True)
    parser.add_argument("--output_csv", type=str, required=True)
    parser.add_argument("--dataset_name", type=str, required=True)
    parser.add_argument("--splitter_name", type=str, required=True)
    parser.add_argument("--featurizer_name", type=str, required=True)
    parser.add_argument("--model_name", type=str, required=True)
    parser.add_argument("--model_path", type=str, required=True)
    args = parser.parse_args()
    main(args)