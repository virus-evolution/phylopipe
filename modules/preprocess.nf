#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

project_dir = projectDir
publish_dev = file(params.publish_dev)


process mask_alignment {
    /**
    * Applies a mask to aligned FASTA
    * @input alignment
    * @output alignment_updated
    * @params mask_file
    */

    input:
    path alignment

    output:
    path "${alignment.baseName}.masked.fa"

    script:
    """
    $project_dir/../bin/add_mask.py \
      --in-alignment ${alignment} \
      --out-alignment "${alignment.baseName}.masked.fa" \
      --mask ${mask_file} \
    """
}


process filter_uk {
    /**
    * Filters UK sequences by metadata to exclude biosample duplicates
    * @input fasta, metadata
    */

    input:
    path fasta
    path metadata

    output:
    path "${fasta.baseName}.filtered.fasta", emit: fasta
    path "${metadata.baseName}.filtered.csv", emit: metadata

    script:
    """
    $project_dir/../bin/filter.py \
            --in-fasta ${fasta} \
            --in-metadata ${metadata} \
            --outgroups ${lineage_splits} \
            --out-fasta "${fasta.baseName}.filtered.fasta" \
            --out-metadata "${metadata.baseName}.filtered.csv" \
            --exclude_true duplicate
    """
}


process clean_fasta_headers {
    /**
    * Cleans up strings in FASTA
    * @input fastas
    */

    input:
    path fasta

    output:
    path "${fasta.baseName}.clean.fa", emit: fasta
    path "${fasta.baseName}.map.csv", emit: map

    script:
    """
    $project_dir/../bin/remove_dodgy_symbols.py \
          --in-fasta ${fasta} \
          --out-fasta "${fasta.baseName}.clean.fa" \
          --out-metadata "${fasta.baseName}.map.csv"
    """
}


process clean_fasta_headers_with_tree {
    /**
    * Cleans up strings in FASTA
    * @input fastas
    */

    input:
    path fasta
    path tree

    output:
    path "${fasta.baseName}.clean.fa", emit: fasta
    path "${fasta.baseName}.map.csv", emit: map
    path "${tree.baseName}.clean.tree", emit: tree

    script:
    """
    $project_dir/../bin/remove_dodgy_symbols.py \
          --in-fasta ${fasta} \
          --in-tree ${tree} \
          --out-fasta "${fasta.baseName}.clean.fa" \
          --out-metadata "${fasta.baseName}.map.csv" \
          --out-tree "${tree.baseName}.clean.tree"
    """
}


process clean_metadata {
    /**
    * Applies cleaned map to metadata
    * @input metadata, map
    */

    input:
    path metadata
    path map

    output:
    path "${metadata.baseName}.clean.csv"

    script:
    """
    $project_dir/../bin/apply_map.py \
          --in-metadata ${metadata} \
          --in-map "${map}" \
          --to-clean ${params.annotations} \
          --out-metadata "${metadata.baseName}.clean.csv"
    """
}


process get_keep_tips {
    /**
    * Pulls out list of tips in metadata that match tree
    * @input metadata, tree
    */

    input:
    path metadata
    path tree

    output:
    path "${metadata.baseName}.tips.txt"

    script:
    """
    fastafunk get_tips \
      --in-metadata ${metadata} \
      --in-tree ${tree} \
      --out-tips "${metadata.baseName}.tips.txt"
    """
}

process prune_tree_with_metadata {
    /**
    * Removes tips of tree with no matching metadata
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
      -r \
      -o "${tree.baseName}.pruned.newick"
    """
    else
    """
    cp ${tree} "${tree.baseName}.pruned.newick"
    """
}

process get_tree_tips {
    /**
    * Gets list of tree tips
    * @input tree
    */

    input:
    path tree

    output:
    path "${tree.baseName}.tips.txt"

    script:
    """
    gotree stats tips \
      -i ${tree} | cut -f4 > "${tree.baseName}.tips.txt"
    """
}


process annotate_metadata {
    /**
    * Adds note column with info about sequences in input
    * @input metadata
    * @output metadata
    */

    input:
    path metadata
    path tips

    output:
    path "${metadata.baseName}.annotated.csv"

    script:
    """
    #!/usr/bin/env python3
    import csv

    tips = set()
    with open("${tips}", 'r', newline = '') as tips_in:
        for line in tips_in:
            tip = line.rstrip().replace("'","")
            tips.add(tip)

    print("Tree contains %d tips" %len(tips))

    with open("${metadata}", 'r', newline = '') as csv_in, \
        open("${metadata.baseName}.annotated.csv", 'w', newline = '') as csv_out:
        reader = csv.DictReader(csv_in, delimiter=",", quotechar='\"', dialect = "unix")
        new_fieldnames = reader.fieldnames
        if "note" not in reader.fieldnames:
            new_fieldnames.append("note")
        writer = csv.DictWriter(csv_out, fieldnames = new_fieldnames, delimiter=",", quotechar='\"', quoting=csv.QUOTE_MINIMAL, dialect = "unix")
        writer.writeheader()
        for row in reader:
            note = []
            if row["note"]:
                note.extend(row["note"].split("|"))
            statement = "tip in input tree"
            if row["sequence_name"] and row["sequence_name"] in tips and statement not in note:
                note.append(statement)
            row["note"] = "|".join(note)
            writer.writerow(row)
    """
}


process announce_summary {
    /**
    * Summarizes subsampling into JSON
    * @input fasta
    */

    input:
    path original
    path deduplicated

    output:
    path "announce.json"

    script:
        if (params.webhook)
            """
            echo '{"text":"' > announce.json
                echo "*${params.whoami}: Preprocessing ${params.date} for tree*\\n" >> announce.json
                echo "> Number of sequences in COG and GISAID input files : \$(cat ${original} | grep '>' | wc -l)\\n" >> announce.json
                echo "> Number of sequences after filtering uk sequences (deduplication by biosample id): \$(cat ${deduplicated} | grep '>' | wc -l)\\n" >> announce.json
                echo '"}' >> announce.json

            echo 'webhook ${params.webhook}'

            curl -X POST -H "Content-type: application/json" -d @announce.json ${params.webhook}
            """
        else
            """
            touch "announce.json"
            """
}

process announce_metadata_pruned_tree {
    /**
    * Announces metadata pruned tree
    * @input tree
    */

    input:
    path tree
    path metadata
    path pruned_tree

    output:
    path "metadata_pruned_tree.json"

    script:
        if (params.webhook)
            """
            echo '{"text":"' > metadata_pruned_tree.json
            echo "*${params.whoami}: Metadata pruned tree for ${params.date} complete*\\n" >> metadata_pruned_tree.json
            echo "> Total number of sequences in original tree: \$(gotree stats tips -i ${tree} | tail -n+2 | wc -l)\\n" >> metadata_pruned_tree.json
            echo "> Total number of sequences in pruned tree: \$(gotree stats tips -i ${pruned_tree} | tail -n+2 | wc -l)\\n" >> metadata_pruned_tree.json
            echo "> Total number of sequences in metadata: \$(tail -n+1 ${metadata} | wc -l)\\n" >> metadata_pruned_tree.json
            echo '"}' >> metadata_pruned_tree.json

            echo 'webhook ${params.webhook}'

            curl -X POST -H "Content-type: application/json" -d @metadata_pruned_tree.json ${params.webhook}
            """
        else
           """
           echo '{"text":"' > metadata_pruned_tree.json
           echo "*${params.whoami}: Metadata pruned tree for ${params.date} complete*\\n" >> metadata_pruned_tree.json
           echo "> Total number of sequences in original tree: \$(gotree stats tips -i ${tree} | tail -n+2 | wc -l)\\n" >> metadata_pruned_tree.json
           echo "> Total number of sequences in pruned tree: \$(gotree stats tips -i ${pruned_tree} | tail -n+2 | wc -l)\\n" >> metadata_pruned_tree.json
           echo "> Total number of sequences in metadata: \$(tail -n+1 ${metadata} | wc -l)\\n" >> metadata_pruned_tree.json
           echo '"}' >> metadata_pruned_tree.json
           """
}

lineage_splits = file(params.lineage_splits)
mask_file = file(params.mask)

workflow mask_and_filter {
    take:
        fasta
        metadata
    main:
        mask_alignment(fasta)
        filter_uk(mask_alignment.out, metadata)
        announce_summary(fasta, filter_uk.out.fasta)
    emit:
        fasta = filter_uk.out.fasta
        metadata = filter_uk.out.metadata
}


workflow clean_fasta_and_metadata {
    take:
        fasta
        metadata
    main:
        clean_fasta_headers(fasta)
        clean_metadata(metadata, clean_fasta_headers.out.map)
    emit:
        fasta = clean_fasta_headers.out.fasta
        metadata = clean_metadata.out
}


workflow clean_fasta_and_metadata_and_tree {
    take:
        fasta
        metadata
        tree
    main:
        clean_fasta_headers_with_tree(fasta, tree)
        clean_metadata(metadata, clean_fasta_headers_with_tree.out.map)
        get_keep_tips(clean_metadata.out, clean_fasta_headers_with_tree.out.tree)
        prune_tree_with_metadata(clean_fasta_headers_with_tree.out.tree, get_keep_tips.out)
        get_tree_tips(prune_tree_with_metadata.out)
        annotate_metadata(clean_metadata.out, get_tree_tips.out)
        announce_metadata_pruned_tree(tree, metadata, prune_tree_with_metadata.out)
    emit:
        fasta = clean_fasta_headers_with_tree.out.fasta
        metadata = annotate_metadata.out
        tree = prune_tree_with_metadata.out
}

workflow {
    fasta = file(params.fasta)
    metadata = file(params.metadata)

    mask_and_filter(fasta, metadata)

    if ( params.newick_tree ) {
        tree = file(params.tree)
        clean_fasta_and_metadata_and_tree(mask_and_filter.out.fasta, mask_and_filter.out.metadata, tree)
    } else {
        clean_fasta_and_metadata(mask_and_filter.out.fasta, mask_and_filter.out.metadata)
    }
}
