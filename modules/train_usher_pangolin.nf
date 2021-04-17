#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

project_dir = projectDir
publish_dir = file(params.publish_dir)
publish_dev = file(params.publish_dev)

//for line in $(cat pangoLEARN/pangoLEARN/data/lineages.downsample.csv); do echo -e "$(echo $line | cut -f2 --delim ',')\t$(echo $line | cut -f1 --delim ',')"; done > lineages.tsv
process make_lineage_annotated_tree {
    /**
    * Publishes master metadata csv for this category
    * @input protobuf, lineages
    * @output protobuf
    */

    publishDir "${publish}/pangolin", pattern: "*.pb", mode: 'copy'

    input:
    path protobuf
    path lineages

    output:
    path "${protobuf.baseName}.lineages.pb"

    script:
    """
    matUtils annotate -i ${protobuf} -c ${lineages} -o "${protobuf.baseName}.lineages.pb"
    """
}

process anonymize_protobuf {
    /**
    * Publishes master metadata csv for this category
    * @input protobuf
    * @output protobuf
    */

    publishDir "${publish}/pangolin", pattern: "*.pb", mode: 'copy'

    input:
    path protobuf

    output:
    path "lineageTree.pb"

    script:
    """
    matUtils mask -i ${protobuf} -S -o lineageTree.pb
    """
}


publish_recipes = file(params.lineage_designations)


workflow train_usher_pangolin {
    take:
        fasta
        metadata
        variants
        newick_tree
        nexus_tree
    main:
        publish_master_metadata(metadata,params.category)
        fetch_min_metadata(fasta,metadata)
        split_recipes(publish_recipes)
        recipe_ch = split_recipes.out.flatten()
        fetch_min_metadata.out.fasta.combine(fetch_min_metadata.out.min_metadata)
                                    .combine(metadata)
                                    .combine(variants)
                                    .combine(newick_tree)
                                    .combine(nexus_tree)
                                    .combine(recipe_ch)
                                    .set{ publish_input_ch }
        publish_tree_recipes(publish_input_ch)
        outputs_ch = publish_tree_recipes.out.flag.collect()
        announce_to_webhook(outputs_ch, "${params.whoami}")
        //publish_s3(publish_tree_recipes.out.tree)
    emit:
        published = outputs_ch
}


workflow {
    fasta = Channel.fromPath(params.fasta)
    metadata = Channel.fromPath(params.metadata)
    variants = Channel.fromPath(params.variants)
    newick_tree = Channel.fromPath(params.newick_tree)
    nexus_tree = Channel.fromPath(params.nexus_tree)

    publish_trees(fasta,
                   metadata,
                   variants,
                   newick_tree,
                   nexus_tree)
}
