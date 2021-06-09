#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

project_dir = projectDir
publish_dev = file(params.publish_dev)


process annotate_tree_uk {
    /**
    * Adds metadata annotations to tree
    * @input tree, metadata
    */

    input:
    path tree
    path metadata

    output:
    path "annotated1.tree"

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
        --output "annotated1.tree"
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
    label 'retry_increasing_mem'

    input:
    path tree
    path metadata

    output:
    path "annotated2.tree"

    script:
    """
    clusterfunk annotate_tips \
        --in-metadata ${metadata} \
        --trait-columns uk_lineage \
        --index-column sequence_name \
        --input ${tree} \
        --output "annotated2.tree"
    """
}

process dequote_tree {
    /**
    * Dequotes tree
    * @input tree
    */

    publishDir "${publish_dev}/trees", pattern: "*.tree", mode: 'copy', saveAs: { "cog_global.${params.date}.newick" }, overwrite: true

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
    publishDir "${publish_dev}/trees", pattern: "*.tree", mode: 'copy', saveAs: { "cog_global.${params.date}.nexus" }, overwrite: true
    label 'retry_increasing_mem'

    input:
    path tree
    path metadata

    output:
    path "annotated3.tree"

    script:
    """
    clusterfunk annotate_tips \
        --in-metadata ${metadata} \
        --trait-columns phylotype \
        --index-column sequence_name \
        --input ${tree} \
        --output "annotated3.tree"
    """
}

process announce_annotation_complete {
    /**
    * Announces expanded and annotated tree
    * @input tree
    */

    input:
    path tree

    output:
    path "full_tree.json"

    script:
        if (params.webhook)
            """
            echo '{"text":"' > full_tree.json
            echo "*${params.whoami}: Annotated tree for ${params.date} complete*\\n" >> full_tree.json
            echo '"}' >> full_tree.json

            echo 'webhook ${params.webhook}'

            curl -X POST -H "Content-type: application/json" -d @full_tree.json ${params.webhook}
            """
        else
           """
           touch "full_tree.json"
           """
}


process annotate_tree {
    /**
    * Adds metadata annotations to tree
    * @input tree, metadata
    */
    label 'retry_increasing_mem'

    input:
    path tree
    path metadata

    output:
    path "annotated.tree"

    script:
    """
    jclusterfunk annotate \
      -i ${tree} \
      -c sequence_name \
      -m ${metadata} \
      --tip-attributes ${params.annotations} \
      -o "annotated.tree" \
      --ignore-missing
    """
}


workflow get_cog_uk_phylotypes {
    take:
        tree
        metadata
    main:
        annotate_tree_uk(tree, metadata)
        deltran_ancestral_reconstruction(annotate_tree_uk.out)
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
        dequote_tree(tree)
        cut_out_tree(dequote_tree.out, uk_lineages_ch)
        phylotype_cut_tree(cut_out_tree.out)
        get_uk_phylotypes_csv(phylotype_cut_tree.out)
        get_uk_phylotypes_csv.out.collectFile(keepHeader: true, skip: 1, name: 'uk_phylotypes.csv')
                                 .set{ uk_phylotypes_csv }
        update_metadata_with_phylotypes(update_uk_lineage_metadata.out,uk_phylotypes_csv)
        annotate_tree_phylotype(annotate_tree_uk_lineage.out, update_metadata_with_phylotypes.out)
    emit:
        nexus_tree = annotate_tree_phylotype.out
        metadata = update_metadata_with_phylotypes.out
        newick_tree = dequote_tree.out
}


workflow post_process_tree {
    take:
        tree
        metadata
    main:
        if ( params.cog_uk ) {
           get_cog_uk_phylotypes(tree, metadata)
           nexus_ch = get_cog_uk_phylotypes.out.nexus_tree
           newick_ch = get_cog_uk_phylotypes.out.newick_tree
           metadata_ch = get_cog_uk_phylotypes.out.metadata
        } else {
           annotate_tree(tree, metadata)
           nexus_ch = annotate_tree.out
           newick_ch = tree
           metadata_ch = metadata

        }
        announce_annotation_complete(nexus_ch)
    emit:
        nexus_tree = nexus_ch
        metadata = metadata_ch
        newick_tree = newick_ch
}

workflow {
    tree = file(params.newick_tree)
    metadata = file(params.metadata)

    post_process_tree(tree,metadata)
}