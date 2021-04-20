#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

project_dir = projectDir
publish_dev = file(params.publish_dev)

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
          --in_map "${map}" \
          --out-metadata "${metadata.baseName}.clean.csv"
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
    path "${fasta.baseName}.tips.fasta", emit: fasta
    path "${fasta.baseName}.new.fasta", emit: to_add

    script:
    """
    fastafunk extract \
        --in-fasta ${fasta} \
        --in-tree ${tree} \
        --out-fasta "${fasta.baseName}.tips.fasta" \
        --reject-fasta "${fasta.baseName}.new.fasta" \
        --low-memory
    """
}

process extract_protected_fasta {
    /**
    * Extracts fasta corresponding to lineage designations fasta
    * @input fasta, tree
    */
    label 'retry_increasing_mem'

    input:
    path fasta
    path metadata

    output:
    path "${fasta.baseName}.protected.fasta"

    script:
    """
    fastafunk extract \
        --in-fasta ${fasta} \
        --in-metadata ${metadata} \
        --out-fasta "${fasta.baseName}.protected.fasta" \
        --reject-fasta "${fasta.baseName}.new.fasta" \
        --low-memory
    """
}

process add_reference_to_fasta {
    /**
    * Creates a new fasta with reference first
    * @input fasta
    */

    input:
    path fasta

    output:
    path "${fasta.baseName}.with_reference.fasta"

    script:
    """
    cat ${reference} > "${fasta.baseName}.with_reference.fasta"
    cat ${fasta} >> "${fasta.baseName}.with_reference.fasta"
    """
}

process fasta_to_vcf {
    /**
    * Makes VCF for usher
    * @input fasta
    */
    memory { 30.0.GB + 10.GB * task.attempt }
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries = 2


    input:
    path fasta

    output:
    path "${fasta.baseName}.vcf"

    script:
    """
    faToVcf ${fasta} ${fasta.baseName}.vcf
    """
}

process usher_start_tree {
    /**
    * Makes usher mutation annotated tree
    * @input tree, vcf
    */
    memory { 30.0.GB + 10.GB * task.attempt }
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries = 1
    cpus 8

    input:
    path vcf
    path tree

    output:
    path "trees/${tree.baseName}.USH.tree", emit: tree
    path "trees/${tree.baseName}.pb", emit: protobuf

    script:
    """
    mkdir -p trees
    usher --tree ${tree} \
          --vcf ${vcf} \
          --threads ${task.cpus} \
          --save-mutation-annotated-tree trees/${tree.baseName}.pb \
          --max-uncertainty-per-sample ${params.max_parsimony_placements} \
          --retain-input-branch-lengths \
          --collapse-tree \
          --write-uncondensed-final-tree \
          --outdir trees

    cp trees/uncondensed-final-tree.nh trees/${tree.baseName}.USH.tree
    """
}

process usher_update_tree {
    /**
    * Makes usher mutation annotated tree
    * @input tree, vcf
    */
    publishDir "${publish_dev}/trees", pattern: "trees/*.pb", mode: 'copy', saveAs: { "cog_global.${params.date}.pb" }, overwrite: true
    maxForks 1
    memory { 50.0.GB + 10.GB * task.attempt }
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'ignore' }
    maxRetries = 2
    cpus 8

    input:
    path vcf
    path protobuf

    output:
    path "trees/${protobuf.baseName}.USH.tree", emit: tree
    path "trees/*.pb", emit: protobuf

    script:
    """
    mkdir -p trees
    usher -i ${protobuf} \
          --vcf ${vcf} \
          --threads ${task.cpus} \
          --save-mutation-annotated-tree trees/${protobuf.baseName}.${params.date}.pb \
          --max-uncertainty-per-sample ${params.max_parsimony_placements} \
          --retain-input-branch-lengths \
          --collapse-tree \
          --write-uncondensed-final-tree \
          --outdir trees
    if [ \$? -eq 0 ]; then
        cp trees/${protobuf.baseName}.${params.date}.pb ${protobuf}
        cp trees/uncondensed-final-tree.nh trees/${protobuf.baseName}.USH.tree
    fi
    """
}

process usher_force_update_tree {
    /**
    * Makes usher mutation annotated tree
    * @input protobuf, vcf
    */
    publishDir "${publish_dev}/trees", pattern: "trees/*.pb", mode: 'copy', saveAs: { "cog_global.${params.date}.pb" }, overwrite: true
    maxForks 1
    memory { 50.0.GB + 10.GB * task.attempt }
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'ignore' }
    maxRetries = 2
    cpus 8

    input:
    path vcf
    path protobuf

    output:
    path "trees/${protobuf.baseName}.USH.tree", emit: tree
    path "trees/*.pb", emit: protobuf

    script:
    """
    mkdir -p trees
    usher -i ${protobuf} \
          --vcf ${vcf} \
          --threads ${task.cpus} \
          --save-mutation-annotated-tree trees/${protobuf.baseName}.${params.date}.pb \
          --retain-input-branch-lengths \
          --collapse-tree \
          --write-uncondensed-final-tree \
          --outdir trees
    if [ \$? -eq 0 ]; then
        cp trees/uncondensed-final-tree.nh trees/${protobuf.baseName}.USH.tree
    fi
    """
}

process root_tree {
    /**
    * Roots tree with WH04
    * @input tree
    */

    input:
    path tree

    output:
    path "${tree.baseName}.rooted.tree"

    script:
    """
    jclusterfunk reroot \
        --format newick \
        -i ${tree} \
        -o ${tree.baseName}.rooted.tree \
        --outgroup Wuhan/WH04/2020
    #touch ${tree.baseName}.rooted.tree
    """
}

process announce_tree_complete {
    /**
    * Announces usher updated tree
    * @input tree
    */

    input:
    path tree

    output:
    path "usher_tree.json"

    script:
        if (params.webhook)
            """
            echo '{"text":"' > usher_tree.json
            echo "*${params.whoami}: Usher expanded tree for ${params.date} complete*\\n" >> usher_tree.json
            echo '"}' >> usher_tree.json

            echo 'webhook ${params.webhook}'

            curl -X POST -H "Content-type: application/json" -d @usher_tree.json ${params.webhook}
            """
        else
           """
           touch "usher_tree.json"
           """
}


process announce_protobuf_complete {
    /**
    * Announces usher updated tree
    * @input protobuf
    */

    input:
    path protobuf

    output:
    path "usher_pb.json"

    script:
        if (params.webhook)
            """
            echo '{"text":"' > usher_pb.json
            echo "*${params.whoami}: Usher protobuf for ${params.date} complete*\\n" >> usher_pb.json
            echo '"}' >> usher_pb.json

            echo 'webhook ${params.webhook}'

            curl -X POST -H "Content-type: application/json" -d @usher_pb.json ${params.webhook}
            """
        else
           """
           touch "usher_pb.json"
           """
}

reference = file(params.reference_fasta)

workflow iteratively_update_tree {
    take:
        fasta
        protobuf
    main:
        fasta.splitFasta( by: params.chunk_size, file: true ).set{ fasta_chunks }
        add_reference_to_fasta(fasta_chunks)
        fasta_to_vcf(add_reference_to_fasta.out)

        usher_update_tree(fasta_to_vcf.out, protobuf)
        final_tree = usher_update_tree.out.tree.last()
        final_protobuf = usher_update_tree.out.protobuf.last()
    emit:
        tree = final_tree
        protobuf = final_protobuf
}

workflow force_update_tree {
    take:
        fasta
        protobuf
        tree
    main:
        lineage_designations = file( params.lineage_designations )
        if ( params.lineage_designations ) {
            lineage_designations = file( params.lineage_designations )
            extract_protected_fasta( fasta, lineage_designations )
            add_reference_to_fasta( extract_protected_fasta.out )
            fasta_to_vcf(add_reference_to_fasta.out)
            usher_force_update_tree(fasta_to_vcf.out, protobuf)
            tree_ch = usher_force_update_tree.out.tree
            protobuf_ch = usher_force_update_tree.out.protobuf
        } else {
            tree_ch = tree
            protobuf_ch = protobuf
        }
    emit:
        tree = tree_ch
        protobuf = protobuf_ch
}


workflow usher_expand_tree {
    take:
        fasta
        metadata
        newick_tree
    main:
        clean_fasta_headers_with_tree(fasta, newick_tree)
        apply_map(metadata, clean_fasta_headers_with_tree.out.map)
        dequote_tree(clean_fasta_headers_with_tree.out.tree)
        extract_tips_fasta(clean_fasta_headers_with_tree.out.fasta, dequote_tree.out)
        add_reference_to_fasta(extract_tips_fasta.out.fasta)
        fasta_to_vcf(add_reference_to_fasta.out)
        usher_start_tree(fasta_to_vcf.out,dequote_tree.out)
        iteratively_update_tree(extract_tips_fasta.out.to_add,usher_start_tree.out.protobuf)
        force_update_tree(clean_fasta_headers_with_tree.out.fasta, iteratively_update_tree.out.protobuf, iteratively_update_tree.out.tree)
        root_tree(force_update_tree.out.tree)
        announce_tree_complete(root_tree.out)
    emit:
        fasta = clean_fasta_headers_with_tree.out.fasta
        metadata = apply_map.out
        tree = root_tree.out
        protobuf = force_update_tree.out.protobuf

}

workflow update_full_tree {
    take:
        fasta
        metadata
        newick_tree
        protobuf
    main:
        clean_fasta_headers_with_tree(fasta, newick_tree)
        apply_map(metadata, clean_fasta_headers_with_tree.out.map)
        dequote_tree(clean_fasta_headers_with_tree.out.tree)
        extract_tips_fasta(clean_fasta_headers_with_tree.out.fasta, dequote_tree.out)
        iteratively_update_tree(extract_tips_fasta.out.to_add,protobuf)
        force_update_tree(clean_fasta_headers_with_tree.out.fasta, iteratively_update_tree.out.protobuf, iteratively_update_tree.out.tree)
        root_tree(force_update_tree.out.tree)
        announce_tree_complete(root_tree.out)
    emit:
        fasta = clean_fasta_headers_with_tree.out.fasta
        metadata = apply_map.out
        tree = root_tree.out
        protobuf = force_update_tree.out.protobuf
}

workflow build_protobuf {
    take:
        fasta
        newick_tree
    main:
        clean_fasta_headers_with_tree(fasta, newick_tree)
        dequote_tree(clean_fasta_headers_with_tree.out.tree)
        add_reference_to_fasta(clean_fasta_headers_with_tree.out.fasta)
        fasta_to_vcf(add_reference_to_fasta.out)
        usher_start_tree(fasta_to_vcf.out,dequote_tree.out)
        announce_protobuf_complete(usher_start_tree.out.protobuf)
    emit:
        protobuf = usher_start_tree.out.protobuf
}

workflow {
    fasta = file(params.fasta)
    newick_tree = file(params.newick_tree)

    usher_expand_tree(fasta,newick_tree)
}