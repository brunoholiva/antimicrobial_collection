#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process RUN_FEATURIZER {
    tag "Featurize: ${dataset.baseName} w/ ${featurizer.baseName}"

    input:
    tuple val(dataset), val(splitter), val(seed), path(train_split), path(test_split), val(featurizer)

    output:
    tuple val(dataset), val(splitter), val(seed), val(featurizer), path("train_features.csv"), path("test_features.csv"), emit: features

    script:

    """

    python ${featurizer} \\
        --train_csv ${train_split} \\
        --test_csv ${test_split} \\
        --train_output_csv train_features.csv \\
        --test_output_csv test_features.csv \\
        --smiles_col ${params.smiles_col} \\
        --activity_col ${params.activity_col}

    """
}