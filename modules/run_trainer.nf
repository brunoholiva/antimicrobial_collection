#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

process RUN_TRAINER {
    tag "Train & CV: ${model.baseName}"

    input:
    tuple val(dataset), val(splitter), val(featurizer), path(features), val(model)

    output:
    tuple val(dataset), val(splitter), val(featurizer), val(model), path("${model.baseName}_${featurizer.baseName}_${splitter.baseName}.joblib"), emit: models
    tuple val(dataset), val(splitter), val(featurizer), val(model), path("cv_results.csv"), emit: cv_results

    script:
    """
    python ${model} \
        --input_csv ${features} \
        --activity_col ${params.activity_col} \
        --output_model_path ${model.baseName}_${featurizer.baseName}_${splitter.baseName}.joblib \
        --output_metrics cv_results.csv \
        --random_state ${params.random_state} \
        --n_jobs ${task.cpus}
    """
}