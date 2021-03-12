#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

project_dir = projectDir
publish_dev = file(params.publish_dev)


process expand_hashmap {
    /**
    * Adds back in identical sequences using hashmap
    * @input tree, hashmap
    */

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

process sort_and_collapse {
    /**
    * Runs gotree to collapse tiny branches
    * @input tree
    */

    input:
    path tree

    output:
    path "cog_gisaid_full.tree"

    script:
    """
    gotree rotate sort -i ${tree} -o rotated.tree
    gotree collapse length --length ${params.collapse} -i rotated.tree -o cog_gisaid_full.tree
    """
}

process annotate_tree {
    /**
    * Adds metadata annotations to tree
    * @input tree, metadata
    */

    input:
    path tree
    path metadata

    output:
    path "annotated.tree"

    script:
    """
    clusterfunk annotate_tips \
        --in-metadata ${metadata} \
        --trait-columns country lineage uk_lineage \
        --index-column sequence_name \
        --boolean-for-trait country='UK' country='UK' \
        --boolean-trait-names country_uk country_uk_deltran \
        --in-format newick \
        --out-format nexus \
        --input ${tree} \
        --output "annotated.tree"
    """
}

process deltran_ancestral_reconstruction {
    /**
    * Infers deltrans
    * @input tree
    */

    input:
    path tree

    output:
    path "${tree.baseName}.deltrans.tree"

    script:
    """
    clusterfunk ancestral_reconstruction \
        --traits country_uk_deltran \
        --deltran \
        --ancestral-state False \
        --input ${tree} \
        --output ${tree.baseName}.deltrans.tree
    """
}

process label_deltran_introductions {
    /**
    * Label deltrans introductions
    * @input tree
    */

    input:
    path tree

    output:
    path "${tree.baseName}.del_labelled.tree"

    script:
    """
    clusterfunk label_transitions \
        --trait country_uk_deltran \
        --to True \
        --transition-name del_introduction \
        --transition-prefix del_trans_ \
        --stubborn \
        --input "${tree}" \
        --output "${tree.baseName}.del_labelled.tree"
    """
}

process merge_sibling_del_introduction {
    /**
    * Merge sibling deltrans introductions
    * @input tree
    */

    input:
    path tree

    output:
    path "${tree.baseName}.del_merged.tree"

    script:
    """
    clusterfunk merge_transitions \
        --trait-to-merge del_introduction \
        --merged-trait-name del_lineage \
        --merge-siblings \
        --input ${tree} \
        --output "${tree.baseName}.del_merged.tree"
    """
}

process output_annotations {
    /**
    * Output annotations
    * @input tree
    */

    input:
    path tree

    output:
    path "all_traits.csv"

    script:
    """
    clusterfunk extract_tip_annotations \
        --traits country uk_lineage del_introduction del_lineage \
        --input ${tree} \
        --output "all_traits.csv"
    """
}

process merge_and_create_new_uk_lineages {
    /**
    * Merge and create new uk lineages
    * @input traits_csv
    */

    input:
    path traits_csv

    output:
    path "updated_traits.csv"

    script:
    """
    $project_dir/../bin/curate_linages.py \
            ${traits_csv} \
            "updated_traits.csv"
    """
}

process generate_sankey_plot {
    /**
    * Make old traits to new traits sankey
    * @input traits_csv, new_traits_csv
    */

    input:
    path traits_csv
    path new_traits_csv

    output:
    path "sankey_links.txt"
    path "sankey.html"

    script:
    """
    $project_dir/../bin/get_sankey_links.py \
            ${traits_csv} \
            ${new_traits_csv} \
            "sankey_links.txt"

    Rscript $project_dir/../bin/plot_sankey.R "sankey_links.txt" "sankey.html"
    """
}

process update_uk_lineage_metadata {
    /**
    * Update metadata with del_lineage del_introduction uk_lineage microreact_lineage
    * @input metadata, traits_csv, uk_lineage_csv
    */

    input:
    path metadata
    path traits_csv
    path uk_lineage_csv

    output:
    path "cog_gisaid.lineages.with_all_traits.csv"

    script:
    """
    fastafunk add_columns \
              --in-metadata ${metadata} \
              --in-data ${traits_csv} \
              --index-column sequence_name \
              --join-on taxon \
              --new-columns del_lineage del_introduction \
              --out-metadata tmp.csv

    fastafunk add_columns \
              --in-metadata tmp.csv \
              --in-data ${uk_lineage_csv} \
              --index-column sequence_name \
              --join-on taxon \
              --new-columns uk_lineage microreact_lineage \
              --out-metadata "cog_gisaid.lineages.with_all_traits.csv"
    """
}

process annotate_tree_uk_lineage {
    /**
    * Adds metadata annotations to tree
    * @input tree, metadata
    */

    input:
    path tree
    path metadata

    output:
    path "annotated.tree"

    script:
    """
    clusterfunk annotate_tips \
        --in-metadata ${metadata} \
        --trait-columns uk_lineage \
        --index-column sequence_name \
        --input ${tree} \
        --output "annotated.tree"
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

process cut_out_tree {
    /**
    * Cuts out tree corresponding to a UK lineage
    * @input tree, tip_file
    */

    input:
    path tree
    path tip_file

    output:
    path "trees/${tip_file.baseName}.tree"

    script:
    """
    mkdir -p trees
    gotree prune -r \
        -i ${tree} \
        --tipfile ${tip_file} \
        -o "trees/${tip_file.baseName}.tree"
    """
}

process phylotype_cut_tree {
    /**
    * Phylotypes the cut out tree corresponding to a UK lineage
    * @input tree
    */

    input:
    path tree

    output:
    path "phylotyped_trees/${tree.baseName}.tree"

    script:
    """
    mkdir -p phylotyped_trees
    clusterfunk phylotype \
        --threshold ${params.phylotype_threshold} \
        --prefix "${tree.baseName}_1" \
        --input ${tree} \
        --in-format newick \
        --output "phylotyped_trees/${tree.baseName}.tree"
    """
}

process get_uk_phylotypes_csv {
    /**
    * Extracts phylotyped csv
    * @input tree
    */

    input:
    path tree

    output:
    path "phylotype_csvs/${tree.baseName}.csv"

    script:
    """
    mkdir -p phylotype_csvs
    clusterfunk extract_tip_annotations \
        --traits country phylotype \
        --input ${tree} \
        --output "phylotype_csvs/${tree.baseName}.csv"
    """
}

process update_metadata_with_phylotypes {
    /**
    * Update metadata with del_lineage del_introduction uk_lineage microreact_lineage
    * @input metadata, traits_csv, uk_lineage_csv
    */

    input:
    path metadata
    path phylotypes_csv

    output:
    path "${metadata.baseName}.with_phylotype_traits.csv"

    script:
    """
    fastafunk add_columns \
              --in-metadata ${metadata} \
              --in-data ${phylotypes_csv} \
              --index-column sequence_name \
              --join-on taxon \
              --new-columns phylotype \
              --out-metadata "${metadata.baseName}.with_phylotype_traits.csv"
    """
}

process annotate_tree_phylotype {
    /**
    * Adds metadata annotations to tree
    * @input tree, metadata
    */

    input:
    path tree
    path metadata

    output:
    path "annotated.tree"

    script:
    """
    clusterfunk annotate_tips \
        --in-metadata ${metadata} \
        --trait-columns phylotype \
        --index-column sequence_name \
        --input ${tree} \
        --output "annotated.tree"
    """
}

workflow post_process_tree {
    take:
        tree
        hashmap
        metadata
    main:
        expand_hashmap(tree,hashmap)
        sort_and_collapse(expand_hashmap.out)
        annotate_tree(sort_and_collapse.out, metadata)
        deltran_ancestral_reconstruction(annotate_tree.out)
        label_deltran_introductions(deltran_ancestral_reconstruction.out)
        merge_sibling_del_introduction(label_deltran_introductions.out)
        output_annotations(merge_sibling_del_introduction.out)

        merge_and_create_new_uk_lineages(output_annotations.out)
        //generate_sankey_plot(output_annotations.out,merge_and_create_new_uk_lineages.out)
        update_uk_lineage_metadata(metadata, output_annotations.out, merge_and_create_new_uk_lineages.out)
        annotate_tree_uk_lineage(merge_sibling_del_introduction.out,update_uk_lineage_metadata.out)

        output_annotations.out.splitCsv(header: true)
                          .map(row -> ["${row.taxon}", "${row.uk_lineage}"])
                          .collectFile() { item -> [ "${item[1]}.txt", "${item[0]}" + '\n' ] }
                          .filter { it.countLines() > 2 }
                          .set{ uk_lineages_ch }
        dequote_tree(sort_and_collapse.out)
        cut_out_tree(dequote_tree.out, uk_lineages_ch)
        phylotype_cut_tree(cut_out_tree.out)
        get_uk_phylotypes_csv(phylotype_cut_tree.out)
        get_uk_phylotypes_csv.out.collectFile(keepHeader: true, skip: 1)
                                 .set{ uk_phylotypes_csv }
        update_metadata_with_phylotypes(update_uk_lineage_metadata.out,uk_phylotypes_csv)
        annotate_tree_phylotype(annotate_tree_uk_lineage.out, update_metadata_with_phylotypes.out)
    emit:
        nexus_tree = annotate_tree_phylotype.out
        metadata = update_metadata_with_phylotypes.out
        newick_tree = dequote_tree.out
        //nexus_tree = annotate_tree_uk_lineage.out
        //metadata = update_uk_lineage_metadata.out
}

workflow {
    tree = file(params.newick_tree)
    hashmap = file(params.hashmap)
    metadata = file(params.metadata)

    post_process_tree(tree,hashmap,metadata)
}