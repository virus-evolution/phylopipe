#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

project_dir = projectDir
publish_dev = file(params.publish_dev)


process usher_tree {
    /**
    * Makes usher mutation annotated tree
    * @input tree, metadata
    */
    publishDir "${publish_dev}", pattern: "trees/*.USH.*", mode: 'copy'

    input:
    path tree
    path fasta

    output:
    path "trees/${tree.baseName}.USH.tree", emit: tree
    path "trees/parsimony-scores.USH.tsv"", emit: scores

    script:
    """
    faToVcf ${fasta} ${fasta.baseName}.vcf
    usher --tree ${tree} \
          --vcf ${fasta.baseName}.vcf \
          --save-mutation-annotated-tree ${tree.baseName}.pb \
          --collapse-tree \
          --write-parsimony-scores-per-node \
          --write-uncondensed-final-tree \
          --outdir trees
    cp trees/uncondensed-final-tree.nh trees/${tree.baseName}.USH.tree
    cp trees/parsimony-scores.tsv trees/parsimony-scores.USH.tsv
    """
}


workflow build_full_tree {
    take:
        fasta
        newick_tree
    main:
        usher_tree(newick_tree,fasta)
    emit:
        tree = usher_tree.out.tree
        scores = usher_tree.out.scores
}

workflow {
    fasta = file(params.fasta)
    newick_tree = file(params.newick_tree)

    build_full_tree(fasta,newick_tree)
}