#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

project_dir = projectDir
publish_dir = file(params.publish_dir)
publish_dev = file(params.publish_dev)


process get_unreliable_tips {
    /**
    * Pulls out list of tips in metadata
    * @input metadata
    */

    input:
    path metadata

    output:
    path "${metadata.baseName}.tips.txt"

    script:
    """
    #!/usr/bin/env python3
    import csv

    filter_column = "is_unreliable_in_tree"
    index_column = "sequence_name"
    with open("${metadata}", 'r', newline = '') as csv_in, \
         open("${metadata.baseName}.tips.txt", 'w') as tips_out:
        reader = csv.DictReader(csv_in, delimiter=",", quotechar='\"', dialect = "unix")
        if index_column not in reader.fieldnames:
            sys.exit("Index column %s not in CSV" %index_column)
        if filter_column in reader.fieldnames:
            for row in reader:
                if row[filter_column] in ["Y","Yes","yes","y",True,"True"]:
                    name = row[index_column].replace('"','').replace("'","")
                    tips_out.write("'%s'\\n" %name)
                    tips_out.write("%s\\n" %name)
    """
}

process prune_unreliable_tips {
    /**
    * Removes unreliable tips of tree
    * @input metadata, tree
    */

    input:
    path tree
    path tips

    output:
    path "${tree.baseName}.pruned.newick"

    script:
    if ( params.prune )
    """
    gotree prune \
      -i ${tree} \
      -f ${tips} \
      -o "${tree.baseName}.pruned.newick"
    """
    else
    """
    cp ${tree} "${tree.baseName}.pruned.newick"
    """
}

process announce_unreliable_pruned_tree {
    /**
    * Announces unreliable pruned tree
    * @input tree, tips, pruned_tree
    */

    input:
    path tree
    path tips
    path pruned_tree

    output:
    path "unreliable_pruned_tree.json"

    script:
        if (params.webhook)
            """
            echo '{"text":"' > unreliable_pruned_tree.json
            echo "*${params.whoami}: Pruned tree with unreliable tips removed for ${params.date} complete*\\n" >> unreliable_pruned_tree.json
            echo "> Total number of sequences in original tree: \$(gotree stats tips -i ${tree} | tail -n+2 | wc -l)\\n" >> unreliable_pruned_tree.json
            echo "> Total number of sequences in pruned tree: \$(gotree stats tips -i ${pruned_tree} | tail -n+2 | wc -l)\\n" >> unreliable_pruned_tree.json
            echo "> Total number of sequences in unreliable list: \$(tail -n+1 ${tips} | wc -l)\\n" >> unreliable_pruned_tree.json
            echo '"}' >> metadata_pruned_tree.json

            echo 'webhook ${params.webhook}'

            curl -X POST -H "Content-type: application/json" -d @unreliable_pruned_tree.json ${params.webhook}
            """
        else
           """
            echo '{"text":"' > unreliable_pruned_tree.json
            echo "*${params.whoami}: Pruned tree with unreliable tips removed for ${params.date} complete*\\n" >> unreliable_pruned_tree.json
            echo "> Total number of sequences in original tree: \$(gotree stats tips -i ${tree} | tail -n+2 | wc -l)\\n" >> unreliable_pruned_tree.json
            echo "> Total number of sequences in pruned tree: \$(gotree stats tips -i ${pruned_tree} | tail -n+2 | wc -l)\\n" >> unreliable_pruned_tree.json
            echo "> Total number of sequences in unreliable list: \$(tail -n+1 ${tips} | wc -l)\\n" >> unreliable_pruned_tree.json
            echo '"}' >> metadata_pruned_tree.json

           """
}

process dequote_tree {
    /**
    * Dequotes tree
    * @input tree
    */

    input:
    path tree

    output:
    path "${tree.baseName}.dequote.tree"

    script:
    """
    sed "s/'//g" ${tree} > "${tree.baseName}.dequote.tree"
    """
}

process extract_tips_fasta {
    /**
    * Extracts fasta corresponding to tips in the tree
    * @input fasta, tree
    */
    label 'retry_increasing_mem'

    input:
    path fasta
    path tree

    output:
    path "${tree.baseName}.tips.fasta", emit: fasta

    script:
    """
    fastafunk extract \
        --in-fasta ${fasta} \
        --in-tree ${tree} \
        --out-fasta "${tree.baseName}.tips.fasta" \
        --reject-fasta "${fasta.baseName}.new.fasta" \
        --low-memory
    """
}

process annotate_metadata {
    /**
    * Adds to note column with info about sequences excluded by date
    * @input metadata
    * @output metadata
    */

    input:
    path metadata
    path full_fasta
    path pruned_fasta

    output:
    path "${metadata.baseName}.annotated.csv"

    script:
    """
    #!/usr/bin/env python3
    import csv

    pruned_seqs = set()
    with open("${pruned_fasta}", 'r', newline = '') as fasta_in:
        for line in fasta_in:
            if line.startswith(">"):
                pruned_seqs.add(line.rstrip()[1:])

    full_seqs = set()
    with open("${full_fasta}", 'r', newline = '') as fasta_in:
        for line in fasta_in:
            if line.startswith(">"):
                if line.rstrip()[1:] not in pruned_seqs:
                    full_seqs.add(line.rstrip()[1:])

    print("Full tree had %d tips and pruned tree had %d extra" %(len(full_seqs), len(pruned_seqs)))

    with open("${metadata}", 'r', newline = '') as csv_in, \
        open("${metadata.baseName}.annotated.csv", 'w', newline = '') as csv_out:
        reader = csv.DictReader(csv_in, delimiter=",", quotechar='\"', dialect = "unix")
        new_fieldnames = reader.fieldnames
        for column in ["note", "why_excluded", "is_excluded"]:
            if column not in reader.fieldnames:
                new_fieldnames.append(column)
        writer = csv.DictWriter(csv_out, fieldnames = new_fieldnames, delimiter=",", quotechar='\"', quoting=csv.QUOTE_MINIMAL, dialect = "unix")
        writer.writeheader()
        for row in reader:
            for column in ["note", "why_excluded", "is_excluded"]:
                if column not in row:
                    row[column] = ""

            if row["sequence_name"] in pruned_seqs:
                row["is_excluded"] = "N"
                if row["why_excluded"] != "":
                    row["note"] += "|" + row["why_excluded"]
                    row["why_excluded"] = ""
            elif row["sequence_name"] in full_seqs:
                row["is_excluded"] = "N"
                if row["why_excluded"] != "":
                    row["note"] += "|" + row["why_excluded"]
                    row["why_excluded"] = "pruned from tree"
            else:
                row["is_excluded"] = "Y"
                if row["why_excluded"] == "":
                    if "tried to add with usher" in row["note"]:
                        row["why_excluded"] = "could not place with usher"
                    elif "sample not recent" in row["note"] or "date_filter" in row and "sample_date older than" in row["date_filter"]:
                        row["why_excluded"] = "sample not recent"
                    else:
                        row["why_excluded"] = row["note"]
            writer.writerow(row)
    """
}

process publish_master_metadata {
    /**
    * Publishes master metadata csv for this category
    * @input metadata
    * @output metadata
    */

    publishDir "${publish_dev}", pattern: "*/*.csv", mode: 'copy', overwrite: true

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
          --filter-column sequence_name is_uk is_cog_uk \
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

    errorStrategy 'retry'
    memory = {8.GB * task.attempt}
    maxRetries = 2
    publishDir "${publish_dir}/", pattern: "**/*.*", mode: 'copy', overwrite: true

    input:
    tuple path(fasta),path(metadata),path(mutations),path(constellations),path(newick_tree),path(pruned_newick_tree),path(nexus_tree),path(recipe)

    output:
    path "${recipe.baseName}.done.txt", emit: flag
    path "public/cog_*_tree.newick", optional: true, emit: tree
    path "**/cog_*.*", emit: all

    script:
    """
    $project_dir/../bin/publish_from_config.py \
      --in-fasta ${fasta} \
      --in-metadata ${metadata} \
      --mutations ${mutations} \
      --constellations ${constellations} \
      --newick-tree ${newick_tree} \
      --pruned-newick-tree ${pruned_newick_tree} \
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

    publishDir "${publish_dev}/", pattern: "s3dir", mode: 'copy'

    input:
    path tree

    output:
    path s3dir

    script:
    """
    mkdir -p s3dir
    cp ${tree} s3dir/cog_global_tree.newick
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

workflow dequote_and_extract_tips_fasta {
    take:
        fasta
        newick_tree
    main:
        dequote_tree(newick_tree)
        extract_tips_fasta(fasta, dequote_tree.out)
    emit:
        fasta = extract_tips_fasta.out
}

workflow dequote_and_extract_tips_fasta2 {
    take:
        fasta
        newick_tree
    main:
        dequote_and_extract_tips_fasta(fasta, newick_tree)
    emit:
        fasta = dequote_and_extract_tips_fasta.out.fasta
}


workflow publish_trees {
    take:
        fasta
        metadata
        mutations
        constellations
        newick_tree
        nexus_tree
    main:
        get_unreliable_tips(metadata)
        prune_unreliable_tips(newick_tree,get_unreliable_tips.out).set{ pruned_newick_tree }
        announce_unreliable_pruned_tree(newick_tree, get_unreliable_tips.out, pruned_newick_tree)
        fetch_min_metadata(fasta,metadata)
        fetch_min_metadata.out.fasta.splitFasta( by: params.chunk_size, file: true ).set{ fasta_chunks }
        dequote_and_extract_tips_fasta(fasta_chunks, pruned_newick_tree)
        dequote_and_extract_tips_fasta.out.collectFile(newLine: false).set{ pruned_tips }
        dequote_and_extract_tips_fasta2(fasta_chunks, newick_tree)
        dequote_and_extract_tips_fasta2.out.collectFile(newLine: false).set{ full_tips }
        annotate_metadata(metadata, full_tips, pruned_tips)
        publish_master_metadata(annotate_metadata.out,params.category)
        split_recipes(publish_recipes)
        recipe_ch = split_recipes.out.flatten()
        fasta.combine(annotate_metadata.out)
                 .combine(mutations)
                 .combine(constellations)
                 .combine(newick_tree)
                 .combine(pruned_newick_tree)
                 .combine(nexus_tree)
                 .combine(recipe_ch)
                 .set{ publish_input_ch }
        publish_tree_recipes(publish_input_ch)
        outputs_ch = publish_tree_recipes.out.flag.collect()
        announce_to_webhook(outputs_ch, "${params.whoami}")
        publish_s3(publish_tree_recipes.out.tree)
    emit:
        published = outputs_ch
}


workflow {
    fasta = Channel.fromPath(params.fasta)
    metadata = Channel.fromPath(params.metadata)
    mutations = Channel.fromPath(params.mutations)
    constellations = Channel.fromPath(params.constellations)
    newick_tree = Channel.fromPath(params.newick_tree)
    nexus_tree = Channel.fromPath(params.nexus_tree)

    publish_trees(fasta,
                   metadata,
                   mutations,
                   constellations,
                   newick_tree,
                   nexus_tree)
}
