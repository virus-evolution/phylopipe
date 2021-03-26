#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// import modules
include { subsample_for_tree } from '../modules/subsample_for_tree.nf'
include { build_split_grafted_tree } from '../modules/build_split_grafted_tree.nf'
include { build_full_tree } from '../modules/build_full_tree.nf'
include { update_full_tree } from '../modules/build_full_tree.nf'
include { post_process_tree } from '../modules/post_process_tree.nf'
include { publish_trees } from '../modules/publish_trees.nf'


workflow {
    ch_fasta = Channel.fromPath(params.fasta)
    ch_metadata = Channel.fromPath(params.metadata)
    ch_variants = Channel.fromPath(params.variants)


    subsample_for_tree(ch_fasta,ch_metadata)
    if ( params.protobuf && params.newick_tree) {
        ch_protobuf = Channel.fromPath(params.protobuf)
        ch_tree = Channel.fromPath(params.newick_tree)
        update_full_tree(subsample_for_tree.out.masked_deduped_fasta, ch_tree, ch_protobuf).set{ ch_full_tree }
    } else if ( params.newick_tree ) {
        ch_tree = Channel.fromPath(params.newick_tree)
        build_full_tree(subsample_for_tree.out.masked_deduped_fasta, ch_tree).set{ ch_full_tree }
    } else {
        build_split_grafted_tree(subsample_for_tree.out.fasta, subsample_for_tree.out.metadata, subsample_for_tree.out.hashmap)
        build_full_tree(subsample_for_tree.out.masked_deduped_fasta, build_split_grafted_tree.out.tree).set{ ch_full_tree }
    }
    post_process_tree(ch_full_tree, subsample_for_tree.out.metadata)
    publish_trees(ch_fasta, post_process_tree.out.metadata, ch_variants, post_process_tree.out.newick_tree, post_process_tree.out.nexus_tree)
}