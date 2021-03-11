#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// import modules
include { subsample_for_tree } from '../modules/subsample_for_tree.nf'
include { build_split_grafted_veryfasttree } from '../modules/build_split_grafted_tree.nf'
include { post_process_tree } from '../modules/post_process_tree.nf'
include { publish_trees } from '../modules/publish_trees.nf'


workflow {
    ch_fasta = Channel.fromPath(params.fasta)
    ch_metadata = Channel.fromPath(params.metadata)
    ch_variants = Channel.fromPath(params.variants)


    subsample_for_tree(ch_fasta,ch_metadata)
    build_split_grafted_veryfasttree(subsample_for_tree.out.fasta, ch_metadata)
    post_process_tree(build_split_grafted_tree.out.tree, subsample_for_tree.out.hashmap, subsample_for_tree.out.metadata)
    publish_trees(ch_fasta, post_process_tree.out.metadata, ch_variants, post_process_tree.out.newick_tree, post_process_tree.out.nexus_tree)
}