#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

project_dir = projectDir
publish_dev = file(params.publish_dev)


process hash_non_unique_seqs {
    /**
    * Subsets a unique set of sequences
    * @input fasta, metadata
    */
    memory { fasta.size() * 2.B  + 2.GB * task.attempt }
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries = 2

    input:
    path fasta

    output:
    path "${fasta.baseName}.unique.fasta", emit: fasta
    path "${fasta.baseName}.hashmap.csv", emit: hashmap

    script:
    """
    $project_dir/../bin/hash_non_unique_seqs.py \
        --in-fasta ${fasta} \
        --out-fasta "${fasta.baseName}.unique.fasta" \
        --out-metadata "${fasta.baseName}.hashmap.csv"
    """
}


process expand_hashmap {
    /**
    * Adds back in identical sequences using hashmap
    * @input tree, hashmap
    */
    publishDir "${publish_dev}/trees", pattern: "*.tree", mode: 'copy'


    input:
    path tree
    path hashmap

    output:
    path "${tree.baseName}.expanded.tree"

    script:
    """
    jclusterfunk insert \
        -i "${tree}" \
        --metadata ${hashmap} \
        --unique-only \
        --ignore-missing \
        --format newick \
        -o "${tree.baseName}.expanded.tree"
    """
}