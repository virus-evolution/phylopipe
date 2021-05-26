#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// import modules
include { start } from '../modules/start.nf'
include { mask_and_filter } from '../modules/preprocess.nf'
include { clean_fasta_and_metadata } from '../modules/preprocess.nf'
include { clean_fasta_and_metadata_and_tree } from '../modules/preprocess.nf'
include { extract_protected_sequences } from '../modules/extract_protected_sequences.nf'
include { subsample_for_tree } from '../modules/subsample_for_tree.nf'
include { build_split_grafted_tree } from '../modules/build_split_grafted_tree.nf'
include { build_protobuf } from '../modules/usher_expand_tree.nf'
include { update_protobuf } from '../modules/usher_expand_tree.nf'
include { train_usher_pangolin } from '../modules/train_usher_pangolin.nf'


workflow {
    start()

    if ( ! params.fasta || ! params.metadata || ! params.lineage_designations ) {
        println("Parameters --fasta, --metadata and --lineage_designations must be provided")
        System.exit(1)
    }
    ch_fasta = Channel.fromPath(params.fasta, checkIfExists: true)
    ch_metadata = Channel.fromPath(params.metadata, checkIfExists: true)
    mask_and_filter(ch_fasta, ch_metadata)

    lineage_designations = file( params.lineage_designations, checkIfExists: true )

    if ( params.newick_tree ) {
        ch_tree = Channel.fromPath(params.newick_tree, checkIfExists: true)
        clean_fasta_and_metadata_and_tree(mask_and_filter.out.fasta, mask_and_filter.out.metadata, ch_tree)
        ch_clean_fasta = clean_fasta_and_metadata_and_tree.out.fasta
        ch_clean_metadata = clean_fasta_and_metadata_and_tree.out.metadata
        ch_clean_tree = clean_fasta_and_metadata_and_tree.out.tree
    } else {
        clean_fasta_and_metadata(mask_and_filter.out.fasta, mask_and_filter.out.metadata)
        ch_clean_fasta = clean_fasta_and_metadata.out.fasta
        ch_clean_metadata = clean_fasta_and_metadata.out.metadata
    }

    if ( params.protobuf ) {
        ch_protobuf = Channel.fromPath(params.protobuf, checkIfExists: true)
    } else {
        if ( ! params.newick_tree ) {
            subsample_for_tree(ch_clean_fasta, ch_clean_metadata)
            build_split_grafted_tree(subsample_for_tree.out.fasta, subsample_for_tree.out.metadata, subsample_for_tree.out.hashmap)
            ch_clean_tree = build_split_grafted_tree.out.tree
        }
        build_protobuf(ch_clean_fasta, ch_clean_tree)
        ch_protobuf = build_protobuf.out.protobuf
    }

    ch_protected = extract_protected_sequences(ch_clean_fasta, ch_clean_metadata)
    update_protobuf(ch_protected, ch_protobuf, ch_clean_metadata)
    train_usher_pangolin(update_protobuf.out.protobuf)
}