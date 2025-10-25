#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process RUN_FEATURIZER {
    tag "Featurize: ${dataset.baseName} w/ ${featurizer.baseName}"

    input:
    tuple val(dataset), val(splitter), path(train_split), path(test_split), val(featurizer)

    output:
    tuple val(dataset), val(splitter), val(featurizer), path("train_features.csv"), path("test_features.csv"), emit: features

    script:
    def featurizer_name = featurizer.baseName

    """
    #!/bin/bash
    
    # Check if this is the stateful featurizer
    if [ "${featurizer_name}" == "desc2D" ]; then
        
        # Run the script ONCE
        python ${featurizer} \\
            --train_csv ${train_split} \\
            --test_csv ${test_split} \\
            --train_output_csv train_features.csv \\
            --test_output_csv test_features.csv \\
            --smiles_col ${params.smiles_col}

    else
        
        # Run stateless featurizers (like morgan.py) TWICE
        python ${featurizer} \\
            --input_csv ${train_split} \\
            --output_csv train_features.csv \\
            --smiles_col ${params.smiles_col}
        
        python ${featurizer} \\
            --input_csv ${test_split} \\
            --output_csv test_features.csv \\
            --smiles_col ${params.smiles_col}
    fi
    """
}