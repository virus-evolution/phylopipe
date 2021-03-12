#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

project_dir = projectDir
publish_dir = file(params.publish_dir)
publish_dev = file(params.publish_dev)


process publish_master_metadata {
    /**
    * Publishes master metadata csv for this category
    * @input metadata
    * @output metadata
    */

    publishDir "${publish_dev}", pattern: "*/*.csv", mode: 'copy'

    input:
    path metadata
    val category

    output:
    path "${category}/${category}_master.csv"

    script:
    """
    mkdir -p ${category}
    cp ${metadata} ${category}/${category}_master.csv
    """
}


process fetch_min_metadata {
    /**
    * Gets `sequence_name` column from METADATA and restricts FASTA to non-omit rows
    * @input fasta, metadata
    * @output fasta, min_metadata
    */

    input:
    path fasta
    path metadata

    output:
    path "cog_gisaid.fa", emit: fasta
    path "cog_gisaid_min.csv", emit: min_metadata

    script:
    """
        fastafunk fetch \
          --in-fasta ${fasta} \
          --in-metadata ${metadata} \
          --index-column sequence_name \
          --filter-column sequence_name
          --out-fasta "cog_gisaid.fa" \
          --out-metadata "cog_gisaid_min.csv" \
          --restrict --low-memory
    """
}


process split_recipes {
    input:
    path recipes

    output:
    path "*.json"

    script:
    """
    #!/usr/bin/env python3
    import json
    i = 0

    with open("${recipes}", 'r') as f:
        recipes = json.load(f)

        for d in recipes:
            for entry in recipes[d]:
                new_recipes = {d:[entry]}
                with open("%i.json" %i, 'w') as handle:
                    json.dump(new_recipes,handle)
                i += 1
    """
}


process publish_tree_recipes {
    /**
    * Publishes subsets of combined FASTA, METADATA and TREES for COG-UK and GISAID
    * @input fasta, metadata, newick_tree, nexus_tree,
    * @params publish_recipes.json
    * @output many
    */

    publishDir "${publish_dir}/", pattern: "*/*.*", mode: 'copy', overwrite: false

    input:
    tuple path(fasta),path(min_metadata),path(metadata),path(variants),path(newick_tree),path(nexus_tree),path(recipe)

    output:
    path "*/cog_*.*"

    script:
    """
    $project_dir/../bin/publish_from_config.py \
      --in-fasta ${fasta} \
      --min-metadata ${min_metadata} \
      --full-metadata ${metadata} \
      --variants ${variants}
      --newick-tree ${newick_tree} \
      --nexus-tree ${nexus_tree} \
      --seed ${params.seed} \
      --recipes ${recipe} \
      --date ${params.date}
    """
}

process announce_to_webhook {
    input:
    file published_files
    val name

    script:
    if (params.webhook)
        """
        echo '{"text":"' > announce.json
        echo "*${name} Complete*\\n" >> announce.json
        echo "> Dev outputs in : ${publish_dev}\\n" >> announce.json
        echo "> Publishable outputs in : ${publish_dir}\\n" >> announce.json
        echo '"}' >> announce.json
        echo 'webhook ${params.webhook}'

        curl -X POST -H "Content-type: application/json" -d @announce.json ${params.webhook}
        """
    else
        """
        touch "announce.json"
        """
}


publish_recipes = file(params.publish_recipes)


workflow publish_trees {
    take:
        fasta
        metadata
        variants
        newick_tree
        nexus_tree
    main:
        publish_master_metadata(metadata,"cog_gisaid")
        fetch_min_metadata(fasta,metadata)
        split_recipes(publish_recipes)
        recipe_ch = split_recipes.out.flatten()
        fetch_min_metadata.out.fasta.combine(fetch_min_metadata.out.min_metadata)
                                    .combine(metadata)
                                    .combine(variants)
                                    .combine(newick_tree)
                                    .combine(nexus_tree)
                                    .set{ publish_input_ch }
        publish_recipes(publish_input_ch)
        outputs_ch = publish_recipes.out.all.collect()
        announce_to_webhook(outputs_ch, "Phylopipe2.0")
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
