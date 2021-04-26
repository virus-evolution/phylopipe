#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// import modules
include { mask_and_filter_and_hash } from '../modules/subsample_for_tree.nf'
include { build_split_grafted_tree } from '../modules/build_split_grafted_tree.nf'
include { build_protobuf } from '../modules/usher_expand_tree.nf'
include { extract_protected_fasta } from '../modules/usher_expand_tree.nf'
include { train_usher_pangolin } from '../modules/train_usher_pangolin.nf'


workflow {

    if ( ! params.fasta || ! params.metadata || ! params.lineage_designations ) {
        println("Parameters --fasta, --metadata and --lineage_designations must be provided")
        System.exit(1)
    }
    ch_fasta = Channel.fromPath(params.fasta, checkIfExists: true)
    ch_metadata = Channel.fromPath(params.metadata, checkIfExists: true)
    lineage_designations = file( params.lineage_designations, checkIfExists: true )
    extract_protected_fasta( ch_fasta, lineage_designations )

    if ( params.protobuf ) {
        ch_protobuf = Channel.fromPath(params.protobuf, checkIfExists: true)
    } else {
        mask_and_filter_and_hash(extract_protected_fasta.out,ch_metadata)
        if ( params.newick_tree ) {
            ch_tree = Channel.fromPath(params.newick_tree, checkIfExists: true)
        } else {
            build_split_grafted_tree(mask_and_filter_and_hash.out.fasta, mask_and_filter_and_hash.out.metadata, mask_and_filter_and_hash.out.hashmap)
            build_split_grafted_tree.out.tree.set{ ch_tree }
        }
        build_protobuf(mask_and_filter_and_hash.out.fasta, ch_tree)
        build_protobuf.out.protobuf.set{ ch_protobuf }
    }

    if ( params.update_protobuf ) {
        update_protobuf(ch_protobuf)
        update_protobuf.out.protobuf.set{ ch_complete_protobuf }
    } else {
        ch_complete_protobuf = ch_protobuf
    }
    train_usher_pangolin(ch_complete_protobuf)
}