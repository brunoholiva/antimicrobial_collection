#!/usr/bin/env nextflow
nextflow.enable.dsl = 2


process RUN_SPLITTER {
    tag "Split ${dataset.baseName} - ${splitter.baseName}"

    input:
    tuple val(dataset), val(splitter), val(seed)

    output:
    tuple val(dataset), val(splitter), val(seed), path("train_split.csv"), path("test_split.csv"), emit: splits

    script:
    """
    python ${splitter} \\
        --input_csv ${dataset} \\
        --train_output_csv train_split.csv \\
        --test_output_csv test_split.csv \\
        --smiles_col ${params.smiles_col} \\
        --activity_col ${params.activity_col} \\
        --random_state ${seed}
    """
}
