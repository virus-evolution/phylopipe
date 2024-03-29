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
include { usher_expand_tree } from '../modules/usher_expand_tree.nf'
include { soft_update_usher_tree } from '../modules/usher_expand_tree.nf'
include { hard_update_usher_tree } from '../modules/usher_expand_tree.nf'
include { post_process_tree } from '../modules/post_process_tree.nf'
include { publish_trees } from '../modules/publish_trees.nf'


workflow {
    start()

    ch_fasta = Channel.fromPath(params.fasta)
    ch_metadata = Channel.fromPath(params.metadata)
    ch_mutations = Channel.fromPath(params.mutations)
    ch_constellations = Channel.fromPath(params.constellations)

    mask_and_filter(ch_fasta,ch_metadata)

    if ( params.newick_tree ) {
        ch_tree = Channel.fromPath(params.newick_tree)
        clean_fasta_and_metadata_and_tree(mask_and_filter.out.fasta, mask_and_filter.out.metadata, ch_tree)
        ch_clean_fasta = clean_fasta_and_metadata_and_tree.out.fasta
        ch_clean_metadata = clean_fasta_and_metadata_and_tree.out.metadata
        ch_clean_tree = clean_fasta_and_metadata_and_tree.out.tree
    } else {
        clean_fasta_and_metadata(mask_and_filter.out.fasta, mask_and_filter.out.metadata)
        ch_clean_fasta = clean_fasta_and_metadata.out.fasta
        ch_clean_metadata = clean_fasta_and_metadata.out.metadata
    }

    extract_protected_sequences(ch_clean_fasta, ch_clean_metadata)
    ch_protected = extract_protected_sequences.out.fasta
    ch_protected_metadata = extract_protected_sequences.out.metadata

    if ( ! params.protobuf ) {
        if ( ! params.newick_tree ) {
            subsample_for_tree(ch_clean_fasta, ch_protected_metadata)
            build_split_grafted_tree(subsample_for_tree.out.fasta, subsample_for_tree.out.metadata, subsample_for_tree.out.hashmap)
            ch_clean_tree = build_split_grafted_tree.out.tree
            ch_fasttree_metadata = subsample_for_tree.out.metadata
        } else {
            ch_fasttree_metadata = ch_protected_metadata
        }

        if ( ! params.skip_usher ){
            usher_expand_tree(ch_clean_fasta, ch_clean_tree, ch_fasttree_metadata)
            ch_protobuf = usher_expand_tree.out.protobuf
            ch_expanded_tree = usher_expand_tree.out.tree
            ch_expanded_metadata = usher_expand_tree.out.metadata
        } else {
            ch_expanded_tree = ch_clean_tree
            ch_expanded_metadata = ch_fasttree_metadata
        }

    } else if ( params.newick_tree && params.update_protobuf ) {
        ch_protobuf_raw = Channel.fromPath(params.protobuf)
        soft_update_usher_tree(ch_clean_fasta, ch_clean_tree, ch_protobuf_raw, ch_protected_metadata)
        ch_protobuf = soft_update_usher_tree.out.protobuf
        ch_expanded_tree = soft_update_usher_tree.out.tree
        ch_expanded_metadata = soft_update_usher_tree.out.metadata

    } else if ( params.newick_tree ) {
        ch_protobuf = Channel.fromPath(params.protobuf)
        ch_expanded_tree = ch_clean_tree
        ch_expanded_metadata = ch_protected_metadata
    }

    if (! params.skip_usher ) {
        hard_update_usher_tree(ch_protected, ch_expanded_tree, ch_protobuf, ch_expanded_metadata)
        ch_full_tree = hard_update_usher_tree.out.tree
        ch_full_metadata = hard_update_usher_tree.out.metadata
    } else {
        ch_full_tree = ch_expanded_tree
        ch_full_metadata = ch_expanded_metadata
    }

    post_process_tree(ch_full_tree, ch_full_metadata)
    publish_trees(ch_clean_fasta, post_process_tree.out.metadata, ch_mutations, ch_constellations, post_process_tree.out.newick_tree, post_process_tree.out.nexus_tree)
}