#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// import modules
include { mask_and_filter_and_hash } from '../modules/subsample_for_tree.nf'
include { build_split_grafted_tree } from '../modules/build_split_grafted_tree.nf'
include { build_protobuf } from '../modules/usher_expand_tree.nf'
include { train_usher_pangolin } from '../modules/train_usher_pangolin.nf'


workflow {
    ch_fasta = Channel.fromPath(params.fasta)
    ch_metadata = Channel.fromPath(params.metadata)
    ch_variants = Channel.fromPath(params.variants)

    if ( params.protobuf ) {
        ch_protobuf = Channel.fromPath(params.protobuf)
    } else {
        mask_and_filter_and_hash(ch_fasta,ch_metadata)
        if ( params.newick_tree ) {
            tree_ch = Channel.fromPath(params.newick_tree)
        } else {
            build_split_grafted_tree(mask_and_filter_and_hash.out.fasta, mask_and_filter_and_hash.out.metadata, mask_and_filter_and_hash.out.hashmap).out.tree.set{ tree_ch }
        }
        build_protobuf(mask_and_filter_and_hash.out.fasta, tree_ch).out.protobuf.set{ protobuf_ch }
    }
    train_usher_pangolin(protobuf_ch)
}