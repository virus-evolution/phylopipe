#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

project_dir = projectDir
publish_dir = file(params.publish_dir)
publish_dev = file(params.publish_dev)


process csv_to_tsv {
    input:
    path csv

    output:
    path "${csv.baseName}.tsv"

    script:
    """
    $project_dir/../bin/lineages_csv_to_tsv.py \
      --in-csv ${csv} \
      --out-tsv "${csv.baseName}.tsv"
    """
}

process make_lineage_annotated_tree {
    /**
    * Publishes master metadata csv for this category
    * @input protobuf, lineages
    * @output protobuf
    */

    publishDir "${publish}/pangolin", pattern: "*.pb", mode: 'copy'

    input:
    path protobuf
    path lineages

    output:
    path "${protobuf.baseName}.lineages.pb"

    script:
    """
    matUtils annotate -i ${protobuf} -c ${lineages} -o "${protobuf.baseName}.lineages.pb"
    """
}

process anonymize_protobuf {
    /**
    * Publishes master metadata csv for this category
    * @input protobuf
    * @output protobuf
    */

    publishDir "${publish}/pangolin", pattern: "*.pb", mode: 'copy'

    input:
    path protobuf

    output:
    path "lineageTree.pb"

    script:
    """
    matUtils mask -i ${protobuf} -S -o lineageTree.pb
    """
}




workflow train_usher_pangolin {
    take:
        protobuf
    main:
        lineage_designations = file(params.lineage_designations , checkIfExists: true)
        csv_to_tsv(lineage_designations)
        make_lineage_annotated_tree(protobuf, csv_to_tsv.out)
        anonymize_protobuf(make_lineage_annotated_tree.out)
    emit:
        protobuf = anonymize_protobuf.out
}


workflow {
    protobuf = Channel.fromPath(params.protobuf)

    train_usher_pangolin(protobuf)
}
