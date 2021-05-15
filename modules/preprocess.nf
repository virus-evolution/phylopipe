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


process prune_tree_with_metadata {
    /**
    * Removes tips of tree with no matching metadata
    * @input metadata, tree
    */

    input:
    path metadata
    path tree

    output:
    path "${tree.baseName}.pruned.newick"

    script:
    if ( params.prune )
    """
    $project_dir/../bin/prune_tree.py \
          --metadata ${metadata} \
          --in-tree "${tree}" \
          --out-tree "${tree.baseName}.pruned.newick"
    """
    else
    """
    cp ${tree} "${tree.baseName}.pruned.newick"
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
        prune_tree_with_metadata(clean_metadata.out, clean_fasta_headers_with_tree.out.tree)
    emit:
        fasta = clean_fasta_headers_with_tree.out.fasta
        metadata = clean_metadata.out
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
