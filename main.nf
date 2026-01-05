#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { RUN_SPLITTER } from "./modules/run_splitter.nf"
include { RUN_FEATURIZER } from "./modules/run_featurizer.nf"
include { RUN_TRAINER } from "./modules/run_trainer.nf"
include { RUN_EVALUATOR } from "./modules/run_evaluator.nf"

workflow {
    ch_datasets = channel.fromPath(params.datasets)
    ch_splitters = channel.fromPath(params.splitters)
    ch_featurizers = channel.fromPath(params.featurizers)
    ch_models = channel.fromPath(params.models)
    ch_metrics = channel.value(file(params.metrics))

    RUN_SPLITTER(ch_datasets.combine(ch_splitters))

    RUN_FEATURIZER(RUN_SPLITTER.out.labeled_data.combine(ch_featurizers))

    RUN_TRAINER(RUN_FEATURIZER.out.features.combine(ch_models))

    ch_eval_input = RUN_TRAINER.out.models.join(RUN_TRAINER.out.cv_results, by: [0,1,2,3])
    
    RUN_EVALUATOR(ch_eval_input, ch_metrics)

    RUN_EVALUATOR.out.clean_results.collectFile(
        name: "final_summary.csv",
        storeDir: "${params.outdir}/summary",
        keepHeader: true,
        sort: true,
        skip: 1 
    )
}