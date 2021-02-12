#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

project_dir = projectDir
publish_dir = file(params.publish_dir)
publish_dev = file(params.publish_dev)


process split_recipes {
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


process publish_recipes {
    /**
    * Publishes subsets of combined FASTA and METADATA for COG-UK and GISAID
    * @input gisaid_unaligned_fasta, gisaid_aligned_fasta, gisaid_trimmed_fasta, combined_fasta,
    * gisaid_metadata, combined_metadata, gisaid_variants, combined_variants
    * @params publish_recipes.json
    * @output many
    */

    publishDir "${publish_dir}/", pattern: "*/*.*", mode: 'copy', overwrite: false

    input:
    tuple path(gisaid_fasta),path(gisaid_metadata),path(gisaid_variants),path(recipe)

    output:
    path "*/gisaid_*.*"

    script:
    """
    $project_dir/../bin/publish_from_config.py \
      --recipes ${recipe} \
      --date ${params.date} \
      --gisaid_fasta ${gisaid_fasta} \
      --gisaid_metadata ${gisaid_metadata} \
      --gisaid_variants ${gisaid_variants}
    """
}


process announce_to_webhook {
    input:
    file published_files

    script:
    if (params.webhook)
        """
        echo '{"text":"' > announce.json
        echo "*Datapipe Complete*\\n" >> announce.json
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

recipes = file(params.publish_recipes)

workflow publish_gisaid {
    take:
        gisaid_fasta
        gisaid_metadata
        gisaid_variants
    main:
        split_recipes()
        recipe_ch = split_recipes.out.flatten()
        gisaid_fasta.combine(gisaid_metadata)
                    .combine(gisaid_variants)
                    .combine(recipe_ch)
                    .set{ publish_input_ch }
        publish_recipes(publish_input_ch)
        outputs_ch = publish_recipes.out.collect()
        announce_to_webhook(outputs_ch)
}


workflow {
    gisaid_fasta = Channel.fromPath(params.gisaid_fasta)
    gisaid_metadata = Channel.fromPath(params.gisaid_metadata)
    gisaid_variants = Channel.fromPath(params.gisaid_variants)

    publish_all(gisaid_fasta,
                gisaid_metadata,
                gisaid_variants)
}
