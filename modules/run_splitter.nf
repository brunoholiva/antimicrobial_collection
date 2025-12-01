#!/usr/bin/env nextflow
nextflow.enable.dsl = 2


process RUN_SPLITTER {
    tag "Fold Assignment: ${dataset.baseName} - ${splitter.baseName}"

    input:
    tuple val(dataset), val(splitter)

    output:
    tuple val(dataset), val(splitter), path("labeled_data.csv"), emit: labeled_data

    script:
    """
    python ${splitter} \
        --input_csv ${dataset} \
        --output_csv labeled_data.csv \
        --smiles_col ${params.smiles_col} \
        --activity_col ${params.activity_col} \
        --random_state ${params.random_state}
    """
}
