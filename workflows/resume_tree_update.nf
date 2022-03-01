#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// import modules
include { extract_protected_sequences } from '../modules/extract_protected_sequences.nf'
include { resume_hard_update_usher_tree } from '../modules/usher_expand_tree.nf'
include { post_process_tree } from '../modules/post_process_tree.nf'
include { publish_trees } from '../modules/publish_trees.nf'


workflow {
    ch_fasta = Channel.fromPath(params.fasta)
    ch_metadata = Channel.fromPath(params.metadata)
    ch_mutations = Channel.fromPath(params.mutations)
    ch_constellations = Channel.fromPath(params.constellations)

    extract_protected_sequences(ch_fasta, ch_metadata)
    ch_protected = extract_protected_sequences.out.fasta
    ch_protected_metadata = extract_protected_sequences.out.metadata

    ch_protobuf = Channel.fromPath(params.protobuf)
    ch_vcfs = Channel.fromPath(params.vcfs)
    ch_vcfs.collect().set{ vcf_list }
    ch_usher_log = Channel.fromPath(params.usher_log)

    resume_hard_update_usher_tree(ch_protected, ch_protobuf, ch_protected_metadata, vcf_list, ch_usher_log)

    ch_full_tree = resume_hard_update_usher_tree.out.tree
    ch_full_metadata = resume_hard_update_usher_tree.out.metadata

    post_process_tree(ch_full_tree, ch_full_metadata)
    publish_trees(ch_fasta, post_process_tree.out.metadata, ch_mutations, ch_constellations, post_process_tree.out.newick_tree, post_process_tree.out.nexus_tree)
}