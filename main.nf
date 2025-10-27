#!/usr/bin/env nextflow
nextflow.enable.dsl = 2


// MODULAR CHEMINFORMATICS PIPELINE
// This Nextflow pipeline orchestrates various cheminformatics tasks such as data splitting, featurization, model training, and evaluation.
// Each task is encapsulated in a separate process, allowing for modularity and reusability of components.

include { RUN_SPLITTER } from "./modules/run_splitter.nf"
include { RUN_FEATURIZER } from "./modules/run_featurizer.nf"
include { RUN_TRAINER } from "./modules/run_trainer.nf"
include { RUN_EVALUATOR } from "./modules/run_evaluator.nf"


workflow {
    ch_datasets = channel.fromPath(params.datasets)
    ch_splitters = channel.fromPath(params.splitters)
    ch_featurizers = channel.fromPath(params.featurizers)
    ch_models = channel.fromPath(params.models)
    ch_evaluators = channel.fromPath(params.evaluators)


    ch_for_splitting = ch_datasets.combine(ch_splitters)
    RUN_SPLITTER(ch_for_splitting)

    ch_for_featurization = RUN_SPLITTER.out.splits.combine(ch_featurizers)
    RUN_FEATURIZER(ch_for_featurization)

    ch_for_training = RUN_FEATURIZER.out.features.combine(ch_models)
    RUN_TRAINER(ch_for_training)

    ch_for_evaluation = RUN_TRAINER.out.models.combine(ch_evaluators)
    RUN_EVALUATOR(ch_for_evaluation)

    RUN_EVALUATOR.out.results_csv.collectFile(
        name: "all_results_summary.csv",
        storeDir: "${params.outdir}/summary",
        skip: 1,
        seed: "dataset,split,featurizer,model,evaluator,seed,test_auc,test_accuracy,test_precision,test_recall,test_average_precision\n",
    )
}
