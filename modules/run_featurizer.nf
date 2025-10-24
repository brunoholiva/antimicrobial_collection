#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

process RUN_FEATURIZER {
    tag "Featurize ${dataset.baseName} - ${splitter.baseName} - ${featurizer.baseName}"

    input:
    tuple val(dataset), val(splitter), path(train_split), path(test_split), val(featurizer)

    output:
    tuple val(dataset), val(splitter), val(featurizer), path("train_features.csv"), path("test_features.csv"), emit: features

    script:
    """
    python ${featurizer} \\
        --input_csv ${train_split} \\
        --output_csv train_features.csv \\
        --smiles_col ${params.smiles_col} \\
    
    python ${featurizer} \\
        --input_csv ${test_split} \\
        --output_csv test_features.csv \\
        --smiles_col ${params.smiles_col} \\
    """
}
