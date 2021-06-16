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

process remove_reference {
    /**
    * Removes reference sequences from fasta
    * @input alignment
    * @output alignment_updated
    * @params reference
    */

    input:
    path fasta

    output:
    path "${fasta.baseName}.noref.fa"

    script:
    """
    $project_dir/../bin/remove_ref.py \
      --in-alignment ${fasta} \
      --out-alignment "${fasta.baseName}.noref.fa" \
      --reference ${reference}
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


process optimize_start_tree {
    /**
    * Runs matOptimize to improve tree
    * @input tree
    */

    publishDir "${publish_dev}/trees", pattern: "*.pb", mode: 'copy', saveAs: { "cog_global.${params.date}.pb" }, overwrite: true

    input:
    path protobuf

    output:
    path "${protobuf.baseName}.optimized.pb"

    script:
    """
    matOptimize -i ${protobuf} \
            -o "${protobuf.baseName}.optimized.pb" \
            -r 100 \
            -T ${params.max_cpus} \
            -s 259200
    """
}

process protobuf_to_tree {
    /**
    * Runs matUtils to convert to tree
    * @input tree
    */

    publishDir "${publish_dev}/trees", pattern: "*.tree", mode: 'copy', saveAs: { "cog_global.${params.date}.tree" }, overwrite: true

    input:
    path protobuf

    output:
    path "${protobuf.baseName}.tree"

    script:
    """
    matUtils extract \
        -i ${protobuf} \
        -t "${protobuf.baseName}.tree"
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

process prune_tree_of_long_branches {
    /**
    * Removes tips of tree corresponding to long branches
    * @input tree
    */

    publishDir "${publish_dev}/trees", pattern: "*.pb", mode: 'copy', saveAs: { "cog_global.${params.date}.pb" }, overwrite: true
    publishDir "${publish_dev}/trees", pattern: "*.tree", mode: 'copy', saveAs: { "cog_global.${params.date}.tree" }, overwrite: true

    input:
    path protobuf

    output:
    path "${protobuf.baseName}.pruned.tree", emit: tree
    path "${protobuf.baseName}.pruned.pb", emit: protobuf

    script:
    """
    matUtils extract -i ${protobuf} \
            --max-parsimony ${params.max_parsimony} \
            --max-branch-length ${params.max_branch_length} \
            --write-tree "${protobuf.baseName}.pruned.tree" \
            --collapse-tree \
            --write-mat "${protobuf.baseName}.pruned.pb"
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

process annotate_metadata {
    /**
    * Adds to note column with info about sequences added
    * @input metadata
    * @output metadata
    */

    input:
    path metadata
    path fasta

    output:
    path "${metadata.baseName}.annotated.csv"

    script:
    """
    #!/usr/bin/env python3
    import csv

    seqs = set()
    with open("${fasta}", 'r', newline = '') as fasta_in:
        for line in fasta_in:
            if line.startswith(">"):
                seqs.add(line.rstrip()[1:])

    print("Fasta to add had %d sequences" %len(seqs))

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
            if row["sequence_name"] in seqs:
                note.append("tried to add with usher")
            if row["note"]:
                note.append(row["note"])
            row["note"] = "|".join(note)
            writer.writerow(row)
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


workflow prepare_fasta {
    take:
        fasta
    main:
        remove_reference(fasta)
        mask_fasta(remove_reference.out)
        add_reference_to_fasta(mask_fasta.out, reference)
        fasta_to_vcf(add_reference_to_fasta.out)
    emit:
        vcf = fasta_to_vcf.out
}



workflow iteratively_update_tree {
    take:
        fasta
        protobuf
        metadata
    main:
        fasta.splitFasta( by: params.chunk_size, file: true ).set{ fasta_chunks }
        prepare_fasta(fasta_chunks)
        prepare_fasta.out.collect().set{ vcf_list }
        usher_update_tree(vcf_list, protobuf)
        annotate_metadata(metadata, fasta)
        add_usher_metadata(usher_update_tree.out.usher_log.last(),annotate_metadata.out)
        prune_tree_of_long_branches(usher_update_tree.out.protobuf.last())

        final_tree = prune_tree_of_long_branches.out.tree
        final_protobuf = prune_tree_of_long_branches.out.protobuf
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
        prepare_fasta(fasta_chunks)
        prepare_fasta.out.collect().set{ vcf_list }
        usher_force_update_tree(vcf_list, protobuf)
        annotate_metadata(metadata, fasta)
        add_usher_metadata(usher_force_update_tree.out.usher_log.last(),annotate_metadata.out)
        prune_tree_of_long_branches(usher_force_update_tree.out.protobuf.last())
        final_tree = prune_tree_of_long_branches.out.tree
        final_protobuf = prune_tree_of_long_branches.out.protobuf
    emit:
        tree = final_tree
        protobuf = final_protobuf
        metadata = add_usher_metadata.out
}


workflow usher_expand_tree {
    take:
        fasta
        newick_tree
        metadata
    main:
        dequote_tree(newick_tree)
        extract_tips_fasta(fasta, dequote_tree.out)
        prepare_fasta(extract_tips_fasta.out.fasta)
        usher_start_tree(prepare_fasta.out,dequote_tree.out)
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
        optimize_start_tree(out_pb)
        protobuf_to_tree(optimize_start_tree.out)
        root_tree(protobuf_to_tree.out)
        rescale_branch_lengths(root_tree.out)
        announce_tree_complete(rescale_branch_lengths.out, "initial")
    emit:
        tree = rescale_branch_lengths.out
        protobuf = optimize_start_tree.out
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
        rescale_branch_lengths(root_tree.out)
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
        rescale_branch_lengths(root_tree.out)
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
        prepare_fasta(extract_tips_fasta.out.fasta)
        usher_start_tree(prepare_fasta.out,dequote_tree.out)
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
        iteratively_update_tree(fasta,protobuf,metadata)
        announce_protobuf_complete(iteratively_update_tree.out.protobuf)
    emit:
        protobuf = iteratively_update_tree.out.protobuf
        metadata = iteratively_update_tree.out.metadata
}

workflow {
    fasta = file(params.fasta)
    newick_tree = file(params.newick_tree)
    metadata = file(params.metadata)

    usher_expand_tree(fasta, newick_tree, metadata)
}