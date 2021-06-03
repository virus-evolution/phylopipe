#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

project_dir = projectDir
publish_dev = file(params.publish_dev)


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


process mask_fasta {
    /**
    * Applies a mask to fasta
    * @input alignment
    * @output alignment_updated
    * @params mask_file
    */

    input:
    path fasta

    output:
    path "${fasta.baseName}.masked.fa"

    script:
    """
    $project_dir/../bin/add_mask.py \
      --in-alignment ${fasta} \
      --out-alignment "${fasta.baseName}.masked.fa" \
      --vcf ${mask_file} \
      --filters nanopore_adapter single_src highly_homoplasic
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
    # hack to make sure nexflow doesn't stall if no new sequences
    if [  -z \$(head -n1 "${fasta.baseName}.new.fasta") ]
    then
        head -n10 "${fasta.baseName}.tips.fasta" > "${fasta.baseName}.new.fasta"
    fi
    if [  -z \$(head -n1 "${fasta.baseName}.tips.fasta") ]
    then
        head -n10 "${fasta.baseName}.new.fasta" > "${fasta.baseName}.tips.fasta"
    fi
    """
}


process copy_protobuf {
    /**
    * Copies protobuf
    * @input protobuf
    */

    input:
    path protobuf

    output:
    path "${protobuf.baseName}.in.pb"

    script:
    """
    cp ${protobuf} "${protobuf.baseName}.in.pb"
    """
}


process add_reference_to_fasta {
    /**
    * Creates a new fasta with reference first
    * @input fasta
    */
    maxForks 100
    errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }

    input:
    path fasta
    path ref

    output:
    path "${fasta.baseName}.with_reference.fasta"

    script:
    """
    cat ${ref} > "${fasta.baseName}.with_reference.fasta"
    cat ${fasta} >> "${fasta.baseName}.with_reference.fasta"
    """
}

process fasta_to_vcf {
    /**
    * Makes VCF for usher
    * @input fasta
    */
    maxForks 100
    memory { 10.GB * task.attempt + fasta.size() * 3.B }
    errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }


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
    publishDir "${publish_dev}/trees", pattern: "trees/*.pb", mode: 'copy', saveAs: { "cog_global.${params.date}.pb" }, overwrite: true
    publishDir "${publish_dev}/trees", pattern: "trees/*.tree", mode: 'copy', saveAs: { "cog_global.${params.date}.tree" }, overwrite: true

    memory { 10.GB * task.attempt + vcf.size() * 3.B }
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 1
    cpus {params.max_cpus}

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
    publishDir "${publish_dev}/trees", pattern: "trees/*.tree", mode: 'copy', saveAs: { "cog_global.${params.date}.tree" }, overwrite: true

    memory { 10.GB * task.attempt + vcf_list.size() * 3.B }
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'ignore' }
    maxRetries 2
    cpus {params.max_cpus}

    input:
    path vcf_list
    path protobuf

    output:
    path "trees/${protobuf.baseName}.USH.tree", emit: tree
    path "trees/*.pb", emit: protobuf
    path "usher.log", emit: usher_log

    script:
    """
    mkdir -p trees
    cp ${protobuf} in.pb

    for vcf in ${vcf_list}
    do
      echo "Adding VCF \$vcf to tree\n" >> update_tree.log
      usher -i in.pb \
          --vcf \$vcf \
          --threads ${task.cpus} \
          --save-mutation-annotated-tree out.pb \
          --max-uncertainty-per-sample ${params.max_parsimony_placements} \
          --write-uncondensed-final-tree \
          --outdir trees 2>> usher.log
      echo "Total number of sequences in tree: \$(gotree stats tips -i trees/uncondensed-final-tree.nh | tail -n+2 | wc -l)\n" >> update_tree.log
      if [ \$? -eq 0 ]; then
        mv out.pb in.pb
      fi
    done
    if [ \$? -eq 0 ]; then
        cp in.pb trees/${protobuf.baseName}.${params.date}.pb
        cp trees/uncondensed-final-tree.nh trees/${protobuf.baseName}.USH.tree
        cat usher.log | grep "Number of parsimony-optimal placements"
    fi
    """
}

process usher_force_update_tree {
    /**
    * Makes usher mutation annotated tree
    * @input protobuf, vcf
    */
    publishDir "${publish_dev}/trees", pattern: "trees/*.pb", mode: 'copy', saveAs: { "cog_global.${params.date}.pb" }, overwrite: true
    publishDir "${publish_dev}/trees", pattern: "trees/*.tree", mode: 'copy', saveAs: { "cog_global.${params.date}.tree" }, overwrite: true

    memory { 10.GB * task.attempt + vcf_list.size() * 3.B }
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'ignore' }
    maxRetries 2
    cpus {params.max_cpus}

    input:
    path vcf_list
    path protobuf

    output:
    path "trees/${protobuf.baseName}.USH.tree", emit: tree
    path "trees/*.pb", emit: protobuf
    path "usher.log", emit: usher_log

    script:
    """
    mkdir -p trees
        cp ${protobuf} in.pb

        for vcf in ${vcf_list}
        do
          echo "Adding VCF \$vcf to tree\n" >> update_tree.log
          usher -i in.pb \
              --vcf \$vcf \
              --threads ${task.cpus} \
              --save-mutation-annotated-tree out.pb \
              --write-uncondensed-final-tree \
              --outdir trees 2>> usher.log
          echo "Total number of sequences in tree: \$(gotree stats tips -i trees/uncondensed-final-tree.nh | tail -n+2 | wc -l)\n" >> update_tree.log
          if [ \$? -eq 0 ]; then
            mv out.pb in.pb
          fi
        done
        if [ \$? -eq 0 ]; then
            cp in.pb trees/${protobuf.baseName}.${params.date}.pb
            cp trees/uncondensed-final-tree.nh trees/${protobuf.baseName}.USH.tree
        fi
    """
}

process add_usher_metadata {
    /**
    * Adds metadata from usher log
    * @input log, metadata
    */

    input:
    path usher_log
    path metadata

    output:
    path "${metadata.baseName}.usher.csv"

    script:
    """
    $project_dir/../bin/parse_usher_log.py \
                --in ${usher_log} \
                --out "usher_log.csv"

    fastafunk add_columns \
                  --in-metadata ${metadata} \
                  --in-data "usher_log.csv" \
                  --index-column sequence_name \
                  --join-on sequence_name \
                  --new-columns parsimony_score num_parsimony_optimal_placements is_unreliable_in_tree \
                  --out-metadata "${metadata.baseName}.usher.csv"
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


process cut_long_branches {
    /**
    * Cut branches whose length is greater than or equal to the given length
    * @input tree
    */

    input:
    path tree

    output:
    path "${tree.baseName}.components.tsv"

    script:
    """
    gotree brlen cut \
        -l ${params.max_branch_length} \
        -i ${tree} \
        -o ${tree.baseName}.components.tsv
    """
}

process extract_tips_to_prune {
    /**
    * Parse connected components for individual components
    * @input tree
    */

    input:
    path tsv

    output:
    path "${tsv.baseName}.tips.txt"

    script:
    """
    #!/usr/bin/env python3

    min_component_size = 2

    with open("${tsv}", 'r') as f, \
         open("${tsv.baseName}.tips.txt", 'w') as g:
        for line in f:
            l = line.rstrip().split('\t')
            pos, component_size, tip_list = l

            if int(component_size) < min_component_size:
                tips = tip_list.split()
                for tip in tips:
                    g.write("%s\\n" %tip)
    """
}

process prune_tree_of_long_branches {
    /**
    * Removes tips of tree corresponding to long branches
    * @input tree
    */

    input:
    path tree
    path tips

    output:
    path "${tree.baseName}.pruned.newick"

    script:
    """
    gotree prune \
      -i ${tree} \
      -f ${tips} \
      -o "${tree.baseName}.pruned.newick"
    """
}

process rescale_branch_lengths {
    /**
    * Scales the integer branch lengths to in [0,1]
    * @input tree
    */

    input:
    path tree

    output:
    path "${tree.baseName}.scaled.tree"

    script:
    """
    jclusterfunk scale \
        --format newick \
        -i ${tree} \
        -o ${tree.baseName}.scaled.tree \
        --factor 0.00003344146 \
        --threshold ${params.collapse}
    """
}

process announce_tree_complete {
    /**
    * Announces usher updated tree
    * @input tree
    */

    input:
    path tree
    val label

    output:
    path "usher_tree.json"

    script:
        if (params.webhook)
            """
            echo '{"text":"' > usher_tree.json
            echo "*${params.whoami}: Usher ${label} expanded tree for ${params.date} complete*\\n" >> usher_tree.json
            echo "Total number of sequences in tree: \$(gotree stats tips -i ${tree} | tail -n+2 | wc -l)\n" >> usher_tree.json
            echo '"}' >> usher_tree.json

            echo 'webhook ${params.webhook}'

            curl -X POST -H "Content-type: application/json" -d @usher_tree.json ${params.webhook}
            """
        else
           """
           echo '{"text":"' > usher_tree.json
           echo "*${params.whoami}: Usher ${label} expanded tree for ${params.date} complete*\\n" >> usher_tree.json
           echo "Total number of sequences in tree: \$(gotree stats tips -i ${tree} | tail -n+2 | wc -l)\n" >> usher_tree.json
           echo '"}' >> usher_tree.json
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
           echo '{"text":"' > usher_pb.json
           echo "*${params.whoami}: Usher protobuf for ${params.date} complete*\\n" >> usher_pb.json
           echo '"}' >> usher_pb.json
           """
}

reference = file(params.reference_fasta)
lineage_splits = file(params.lineage_splits)
mask_file = file(params.mask_vcf)


workflow iteratively_update_tree {
    take:
        fasta
        protobuf
        metadata
    main:
        fasta.splitFasta( by: params.chunk_size, file: true ).set{ fasta_chunks }
        add_reference_to_fasta(fasta_chunks, reference)
        mask_fasta(add_reference_to_fasta.out)
        fasta_to_vcf(mask_fasta.out)
        fasta_to_vcf.out.collect().set{ vcf_list }
        usher_update_tree(vcf_list, protobuf)
        add_usher_metadata(usher_update_tree.out.usher_log.last(),metadata)
        final_tree = usher_update_tree.out.tree.last()
        final_protobuf = usher_update_tree.out.protobuf.last()
    emit:
        tree = final_tree
        protobuf = final_protobuf
        metadata = add_usher_metadata.out
}

workflow iteratively_force_update_tree {
    take:
        fasta
        protobuf
        metadata
    main:
        fasta.splitFasta( by: params.chunk_size, file: true ).set{ fasta_chunks }
        add_reference_to_fasta(fasta_chunks, reference)
        mask_fasta(add_reference_to_fasta.out)
        fasta_to_vcf(mask_fasta.out)
        fasta_to_vcf.out.collect().set{ vcf_list }
        usher_force_update_tree(vcf_list, protobuf)
        add_usher_metadata(usher_force_update_tree.out.usher_log.last(),metadata)
        final_tree = usher_force_update_tree.out.tree.last()
        final_protobuf = usher_force_update_tree.out.protobuf.last()
    emit:
        tree = final_tree
        protobuf = final_protobuf
        metadata = add_usher_metadata.out
}

workflow remove_long_branches {
    take:
        tree
    main:
        cut_long_branches(tree)
        extract_tips_to_prune(cut_long_branches.out)
        prune_tree_of_long_branches(tree, extract_tips_to_prune.out)
    emit:
        tree = prune_tree_of_long_branches.out
}


workflow usher_expand_tree {
    take:
        fasta
        newick_tree
        metadata
    main:
        dequote_tree(newick_tree)
        extract_tips_fasta(fasta, dequote_tree.out)
        add_reference_to_fasta(extract_tips_fasta.out.fasta, reference)
        mask_fasta(add_reference_to_fasta.out.fasta)
        fasta_to_vcf(mask_fasta.out)
        usher_start_tree(fasta_to_vcf.out,dequote_tree.out)
        if ( params.update_protobuf ){
            iteratively_update_tree(extract_tips_fasta.out.to_add,usher_start_tree.out.protobuf,metadata)
            out_pb = iteratively_update_tree.out.protobuf
            out_tree = iteratively_update_tree.out.tree
            out_metadata = iteratively_update_tree.out.metadata
        } else {
            out_pb = usher_start_tree.out.protobuf
            out_tree = usher_start_tree.out.tree
            out_metadata = metadata
        }
        root_tree(out_tree)
        remove_long_branches(root_tree.out)
        rescale_branch_lengths(remove_long_branches.out)
        announce_tree_complete(rescale_branch_lengths.out, "initial")
    emit:
        tree = rescale_branch_lengths.out
        protobuf = out_pb
        metadata = out_metadata
}

workflow soft_update_usher_tree {
    take:
        fasta
        newick_tree
        protobuf
        metadata
    main:
        dequote_tree(newick_tree)
        extract_tips_fasta(fasta, dequote_tree.out)
        iteratively_update_tree(extract_tips_fasta.out.to_add,protobuf,metadata)
        root_tree(iteratively_update_tree.out.tree)
        remove_long_branches(root_tree.out)
        rescale_branch_lengths(remove_long_branches.out)
        announce_tree_complete(rescale_branch_lengths.out, "soft")
    emit:
        tree = rescale_branch_lengths.out
        protobuf = iteratively_update_tree.out.protobuf
        metadata = iteratively_update_tree.out.metadata
}

workflow hard_update_usher_tree {
    take:
        fasta
        newick_tree
        protobuf
        metadata
    main:
        dequote_tree(newick_tree)
        extract_tips_fasta(fasta, dequote_tree.out)
        iteratively_force_update_tree(extract_tips_fasta.out.to_add, protobuf,metadata)
        root_tree(iteratively_force_update_tree.out.tree)
        remove_long_branches(root_tree.out)
        rescale_branch_lengths(remove_long_branches.out)
        announce_tree_complete(rescale_branch_lengths.out, "hard")
    emit:
        tree = rescale_branch_lengths.out
        protobuf = iteratively_force_update_tree.out.protobuf
        metadata = iteratively_force_update_tree.out.metadata
}

workflow build_protobuf {
    take:
        fasta
        newick_tree
    main:
        dequote_tree(newick_tree)
        extract_tips_fasta(fasta, dequote_tree.out)
        add_reference_to_fasta(extract_tips_fasta.out.fasta, reference)
        mask_fasta(add_reference_to_fasta.out.fasta)
        fasta_to_vcf(mask_fasta.out)
        usher_start_tree(fasta_to_vcf.out,dequote_tree.out)
        announce_protobuf_complete(usher_start_tree.out.protobuf)
    emit:
        protobuf = usher_start_tree.out.protobuf
}

workflow update_protobuf {
    take:
        fasta
        protobuf
        metadata
    main:
        iteratively_force_update_tree(fasta,protobuf,metadata)
        announce_protobuf_complete(iteratively_force_update_tree.out.protobuf)
    emit:
        protobuf = iteratively_force_update_tree.out.protobuf
        metadata = iteratively_force_update_tree.out.metadata
}

workflow {
    fasta = file(params.fasta)
    newick_tree = file(params.newick_tree)
    metadata = file(params.metadata)

    usher_expand_tree(fasta, newick_tree, metadata)
}