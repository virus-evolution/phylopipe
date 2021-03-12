#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

project_dir = projectDir
publish_dev = file(params.publish_dev)


process usher_tree {
    /**
    * Makes usher mutation annotated tree
    * @input tree, metadata
    */
    label 'retry_increasing_mem'

    publishDir "${publish_dev}", pattern: "trees/*.USH.tsv", mode: 'copy'

    input:
    path tree
    path fasta

    output:
    path "trees/${tree.baseName}.USH.tree", emit: tree
    path "trees/parsimony-scores.USH.tsv", emit: scores

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

process root_tree {
    /**
    * Roots tree with WH04
    * @input tree
    */

    publishDir "${publish_dev}", pattern: "trees/*.tree", mode: 'copy'


    input:
    path tree

    output:
    path "${tree.baseName}.rooted.tree"

    script:
    """
    jclusterfunk reroot \
        --format newick \
        -i ${tree} \
        -o ${tree.baseName}.rooted.tree \
        --outgroup Wuhan/WH04/2020
    #touch ${tree.baseName}.rooted.tree
    """
}

process announce_tree_complete {
    /**
    * Announces usher updated tree
    * @input tree
    */

    input:
    path tree

    output:
    path "usher_tree.json"

    script:
        if (params.webhook)
            """
            echo "{{'text':'" > usher_tree.json
            echo "*Phylopipe2.0: Usher expanded tree for ${params.date} complete*\\n" >> usher_tree.json
            echo "'}}" >> usher_tree.json

            echo 'webhook ${params.webhook}'

            curl -X POST -H "Content-type: application/json" -d @usher_tree.json ${params.webhook}
            """
        else
           """
           touch "usher_tree.json"
           """
}


workflow build_full_tree {
    take:
        fasta
        newick_tree
    main:
        usher_tree(newick_tree,fasta)
        root_tree(usher_tree.out.tree)
        announce_tree_complete(root_tree.out)
    emit:
        tree = root_tree.out.tree
        scores = usher_tree.out.scores
}

workflow {
    fasta = file(params.fasta)
    newick_tree = file(params.newick_tree)

    build_full_tree(fasta,newick_tree)
}