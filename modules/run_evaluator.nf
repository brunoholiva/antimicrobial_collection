process RUN_EVALUATOR {

    tag "Evaluate: ${dataset.baseName} - ${splitter.baseName} - ${featurizer.baseName} - ${model.baseName} - ${evaluator.baseName}"

    publishDir "${params.outdir}/${dataset.baseName}_${splitter.baseName}_${featurizer.baseName}_${model.baseName}", mode: 'copy'

    input:
    tuple val(dataset), val(splitter), val(featurizer), val(model), path(test_features), path(trained_model), val(evaluator)

    output:
    path ("results.csv"), emit: results_csv
    path (trained_model), emit: published_model

    script:

    def dataset_name = dataset.baseName
    def splitter_name = splitter.baseName
    def featurizer_name = featurizer.baseName
    def model_name = model.baseName
    def evaluator_name = evaluator.baseName

    """
    python ${evaluator} \\
        --model_path ${trained_model} \\
        --test_csv ${test_features} \\
        --activity_col ${params.activity_col} \\
        --output_results results.csv \\
        --random_state ${params.random_state} \\
        --dataset_name ${dataset_name} \\
        --split_name ${splitter_name} \\
        --featurizer_name ${featurizer_name} \\
        --model_name ${model_name} \\
        --evaluator_name ${evaluator_name}
    """
}
