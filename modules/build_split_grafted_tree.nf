#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { expand_hashmap } from '../modules/hash.nf'


project_dir = projectDir
publish_dev = file(params.publish_dev)


process get_aliases {
    /**
    * If no alias file provided, pulls latest
    */

    output:
    path "alias_key.json"

    script:
        if (params.lineage_aliases)
            """
            cp ${lineage_aliases} "alias_key.json"
            """
        else
            """
            wget https://raw.githubusercontent.com/cov-lineages/pango-designation/master/pango_designation/alias_key.json
            """
}


process split_fasta {
    /**
    * Splits input fasta for subtrees
    * @input fasta, metadata
    */
    label 'retry_increasing_mem'

    input:
    path fasta
    path metadata
    path aliases

    output:
    path "*.fasta", emit: fasta
    path "*.json", emit: json

    script:
    """
    fastafunk split \
        --in-fasta ${fasta} \
        --in-metadata ${metadata} \
        --index-column sequence_name \
        --index-field lineage \
        --lineage-csv ${lineage_splits} \
        --aliases ${aliases}

    echo '{"text":"' > pre_tree.json
    echo "*${params.whoami}: Ready for ${params.date} tree building*\\n" >> pre_tree.json
    num_lineages=\$(cat ${lineage_splits} | wc -l)
    range={\$num_lineages..1}
    for i in \$(eval echo \${range})
    do
        line=\$(tail -n\$i .command.log | head -n1)
        echo ">\$line\\n" >> pre_tree.json
    done
    echo '"}' >> pre_tree.json
    """
}

process announce_split {
    /**
    * Summarizes splitting
    * @input json
    */

    input:
    path json

    script:
    """
    echo 'webhook ${params.webhook}'
    curl -X POST -H "Content-type: application/json" -d @${json} ${params.webhook}
    """
}

process check_tree_building_size {
    /**
    * Checks fasta size is less than parameter
    * @input fasta
    */

    input:
    path fasta

    output:
    path "${fasta.baseName}.checked.fa"

    script:
    """
    if [[ \$(cat ${fasta} | grep ">" | wc -l) -ge ${params.max_tree_size} ]]
    then
        echo "FASTA TOO BIG TO BUILD"
    else
        cp ${fasta} "${fasta.baseName}.checked.fa"
    fi
    """
}

process fasttree {
    /**
    * Runs fasttree on a lineage
    * @input lineage_fasta
    */
    memory { 4.0 * task.attempt + lineage_fasta.size() * 0.000000065 < 96 ? 4.GB * task.attempt + lineage_fasta.size() * 65.B : 96.GB }
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries = 1
    cpus 3

    input:
    tuple val(lineage), path(lineage_fasta)

    output:
    tuple val(lineage), path("${lineage_fasta.baseName}.unrooted.tree")

    script:
    """
    export OMP_NUM_THREADS=${task.cpus}
    FastTreeMP -nosupport -nt ${lineage_fasta} > ${lineage_fasta.baseName}.unrooted.tree
    #touch ${lineage_fasta.baseName}.unrooted.tree
    """
}

process veryfasttree {
    /**
    * Runs fasttree on a lineage
    * @input lineage_fasta
    */
    memory { 8.0 * task.attempt + lineage_fasta.size() * 0.000000065 < 96 ? 8.GB * task.attempt + lineage_fasta.size() * 65.B : 96.GB }
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries = 1
    cpus 8
    time '2d'

    input:
    tuple val(lineage), path(lineage_fasta)

    output:
    tuple val(lineage), path("${lineage_fasta.baseName}.unrooted.tree")

    script:
    """
    VeryFastTree -double-precision -nosupport -nt ${lineage_fasta} -threads ${task.cpus} > ${lineage_fasta.baseName}.unrooted.tree
    #touch ${lineage_fasta.baseName}.unrooted.tree
    """
}

process root_tree {
    /**
    * Roots tree with given outgroups
    * @input lineage, tree, outgroup
    */

    input:
    tuple val(lineage), path(lineage_tree), val(outgroup)

    output:
    val lineage, emit: lineages
    path "${lineage}.rooted.tree", emit: trees

    script:
    """
    jclusterfunk reroot \
        --format newick \
        -i ${lineage_tree} \
        -o ${lineage}.rooted.tree \
        --outgroup ${outgroup}
    #touch ${lineage}.rooted.tree
    """
}

process graft_tree {
    /**
    * Grafts lineages trees together
    * @input scions, lineages, guide_tree
    */
    memory { 10.GB * task.attempt + scions.size() * 20.B }


    input:
    path scions
    val lineages
    val label

    output:
    path "cog_gisaid_grafted.${label}.tree"

    script:
    """
    echo "${lineages}"
    echo "${scions}"
    clusterfunk graft \
        --scions ${scions} \
        --input ${guide_tree} \
        --output "cog_gisaid_grafted.${label}.tree" \
        --in-format newick \
        --out-format newick \
        --scion-annotation-name scion_lineage \
        --annotate-scions ${lineages}
    #touch cog_gisaid_grafted.tree
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


process announce_tree_complete {
    /**
    * Summarizes splitting
    * @input tree
    */

    input:
    path tree

    output:
    path "grafted_tree.json"

    script:
        if (params.webhook)
            """
            echo '{"text":"' > grafted_tree.json
            echo "*${params.whoami}: Initial grafted tree for ${params.date} complete*\\n" >> grafted_tree.json
            echo "Total number of sequences in tree: \$(gotree stats tips -i ${tree} | tail -n+2 | wc -l)\n" >> grafted_tree.json
            echo '"}' >> grafted_tree.json

            echo 'webhook ${params.webhook}'

            curl -X POST -H "Content-type: application/json" -d @grafted_tree.json ${params.webhook}
            """
        else
            """
            echo '{"text":"' > grafted_tree.json
            echo "*${params.whoami}: Initial grafted tree for ${params.date} complete*\\n" >> grafted_tree.json
            echo "Total number of sequences in tree: \$(gotree stats tips -i ${tree} | tail -n+2 | wc -l)\n" >> grafted_tree.json
            echo '"}' >> grafted_tree.json
            """
}


lineage_splits = file(params.lineage_splits)
guide_tree = file(params.guide_tree)

workflow build_split_grafted_tree {
    take:
        fasta
        metadata
        hashmap
    main:
        if (params.lineage_aliases)
            lineage_aliases = file(params.lineage_aliases)
        get_aliases()
        split_fasta(fasta,metadata,get_aliases.out)
        if (params.webhook)
            announce_split(split_fasta.out.json)
        split_fasta.out.fasta.flatMap { f -> f }.set{ fasta_ch }
        check_tree_building_size(fasta_ch)
        split_fasta.out.fasta.flatMap { f -> f }.map { f -> [f.baseName,f] }.set{ split_fasta_ch }
        Channel.from(lineage_splits).splitCsv(header: false, skip: 1).set{ split_outgroup_ch }

        if ( params.tree_builder == "fasttree" ) {
            fasttree(split_fasta_ch)
            fasttree.out.join(split_outgroup_ch).set{ unrooted_tree_ch }
            label_ch = Channel.from("FT")
        } else if ( params.tree_builder == "veryfasttree" ) {
            veryfasttree(split_fasta_ch)
            veryfasttree.out.join(split_outgroup_ch).set{ unrooted_tree_ch }
            label_ch = Channel.from("VFT")
        } else {
            println "Parameter tree_builder must be either fasttree or veryfasttree, not '${params.tree_builder}'"
            exit 1
        }

        root_tree(unrooted_tree_ch)
        root_tree.out.lineages.toSortedList().set{ lineages_ch }
        root_tree.out.trees.collect().set{ trees_ch }
        graft_tree(trees_ch, lineages_ch, label_ch)
        expand_hashmap(graft_tree.out, hashmap)
        announce_tree_complete(expand_hashmap.out)
    emit:
        tree = expand_hashmap.out
}

workflow finish_split_grafted_tree {
    take:
        fasta
        metadata
        hashmap
    main:
        Channel.from(lineage_splits).splitCsv(header: false, skip: 1).set{ split_outgroup_ch }
        unrooted_trees = Channel.fromPath( "${params.tree_dir}/*.tree" ).map { n -> [ n.getBaseName().replace(".unrooted",""), n ] }
        unrooted_trees.join(split_outgroup_ch).set{ unrooted_tree_ch }
        root_tree(unrooted_tree_ch)
        root_tree.out.lineages.toSortedList().set{ lineages_ch }
        root_tree.out.trees.collect().set{ trees_ch }
        label_ch = Channel.from("FT")
        graft_tree(trees_ch, lineages_ch, label_ch)
        expand_hashmap(graft_tree.out, hashmap)
        sort_and_collapse(expand_hashmap.out)
        announce_tree_complete(sort_and_collapse.out)
    emit:
        tree = sort_and_collapse.out
}


workflow {
    fasta = file(params.unique_fasta)
    metadata = file(params.metadata)
    hashmap = file(params.hashmap)

    if (params.tree_dir) {
        finish_split_grafted_tree(fasta,metadata,hashmap)
    } else {
        build_split_grafted_tree(fasta,metadata,hashmap)
    }
}