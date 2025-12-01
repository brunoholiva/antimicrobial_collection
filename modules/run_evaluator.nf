#!/usr/bin/env nextflow

process RUN_EVALUATOR {
    tag "Parse: ${model.baseName}"

    input:
    tuple val(dataset), val(splitter), val(featurizer), val(model), path("raw_results.csv")
    path parser_script

    output:
    path("clean_results.csv"), emit: clean_results

    script:
    """
    python ${parser_script} \
        --input_csv raw_results.csv \
        --output_csv clean_results.csv \
        --dataset_name ${dataset.baseName} \
        --splitter_name ${splitter.baseName} \
        --featurizer_name ${featurizer.baseName} \
        --model_name ${model.baseName}
    """
}