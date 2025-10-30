#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

process RUN_TRAINER {
    tag "Train: ${dataset.baseName} - ${splitter.baseName} - ${featurizer.baseName} - ${model.baseName}"

    input:
    tuple val(dataset), val(splitter), val(seed), val(featurizer), path(train_features), path(test_features), val(model)

    output:
    tuple val(dataset), val(splitter), val(seed), val(featurizer), val(model), path(test_features), path("trained_model.joblib"), emit: models

    script:
    """
    python ${model} \\
        --train_csv ${train_features} \\
        --activity_col ${params.activity_col} \\
        --output_model trained_model.joblib \\
        --random_state ${seed} \\
        --n_jobs ${task.cpus}
    """
}
