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
    path "cog_gisaid_min.fa", emit: fasta
    path "cog_gisaid_min.csv", emit: min_metadata

    script:
    """
        fastafunk fetch \
          --in-fasta ${fasta} \
          --in-metadata ${metadata} \
          --index-column sequence_name \
          --filter-column sequence_name \
          --out-fasta "cog_gisaid_min.fa" \
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
    path "${recipe.baseName}.done.txt", emit: flag
    path "public/cog_*_tree.newick", optional: true, emit: tree
    path "*/cog_*.*", emit: all

    script:
    """
    $project_dir/../bin/publish_from_config.py \
      --in-fasta ${fasta} \
      --min-metadata ${min_metadata} \
      --full-metadata ${metadata} \
      --variants ${variants} \
      --newick-tree ${newick_tree} \
      --nexus-tree ${nexus_tree} \
      --seed ${params.seed} \
      --recipes ${recipe} \
      --date ${params.date}
    touch "${recipe.baseName}.done.txt"
    """
}

process publish_s3 {
    /**
    * Publishes public files to s3
    * @input tree
    */

    input:
    path tree

    script:
    """
    mkdir -p s3dir
    cp ${tree} s3dir/cog_global_tree.newick

    s3cmd sync s3dir/ s3://cog-uk/phylogenetics/{params.date}/ --acl-public
    s3cmd sync s3dir/ s3://cog-uk/phylogenetics/latest/ --acl-public
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
