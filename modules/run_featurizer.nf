#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process RUN_FEATURIZER {
    tag "Featurize: ${dataset.baseName}"

    input:
    tuple val(dataset), val(splitter), path(input_data), val(featurizer)

    output:
    tuple val(dataset), val(splitter), val(featurizer), path("features.csv"), emit: features

    script:
    """
    python ${featurizer} \
        --input_csv ${input_data} \
        --output_csv features.csv \
        --smiles_col ${params.smiles_col}
    """
}