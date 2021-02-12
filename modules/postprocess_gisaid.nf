#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

project_dir = projectDir
publish_dev = file(params.publish_dev)


process gisaid_counts_by_country {
    /**
    * Outputs number of sequences per country
    * @input metadata
    * @output "gisaid_counts_by_country.csv"
    */
    publishDir "${publish_dev}/${params.category}", pattern: "*.csv", mode: 'copy'


    input:
    path metadata

    output:
    path "gisaid_counts_by_country.csv"

    script:
    """
    fastafunk count \
        --in-metadata ${metadata} \
        --group-column edin_admin_0 \
        --log-file "gisaid_counts_by_country.csv"
    """
}


process gisaid_distance_QC {
    /**
    * Outputs number of sequences per country
    * @input fasta, metadata
    * @output "QC_distances.tsv"
    */
    publishDir "${publish_dev}/${params.category}", pattern: "*.csv", mode: 'copy'


    input:
    path fasta
    path metadata

    output:
    path "QC_distances.tsv"

    script:
    """
    datafunk distance_to_root \
        --input-fasta ${fasta} \
        --input-metadata ${metadata}

    mv distances.tsv "QC_distances.tsv"
    """
}


process gisaid_filter_on_distance_to_WH04 {
    /**
    * Restricts to samples within distance x of WH04
    * @input fasta, metadata, distances
    * @output "QC_distances.tsv"
    */
   publishDir "${publish_dev}/${params.category}", pattern: "*.csv", mode: 'copy', saveAs: {"${params.category}_master.csv"}

    input:
    path fasta
    path metadata
    path distances

    output:
    path "${fasta.baseName}.distance_filtered.fa", emit: fasta
    path "${metadata.baseName}.distance_filtered.csv", emit: metadata

    script:
    """
    #!/usr/bin/env python3
    from Bio import SeqIO
    import csv

    reject = []
    with open("${distances}", 'r', newline = '') as distances_in:
        reader = csv.DictReader(distances_in, delimiter="\t", quotechar='\"', dialect = "unix")
        for row in reader:
            sequence_name = row['sequence_name']
            distance = float(row['distance_stdevs'])
            if distance >= 4.0:
                reject.append(sequence_name)

    alignment = SeqIO.index("${fasta}", "fasta")

    with open("${metadata}", 'r', newline = '') as csv_in, \
        open("${metadata.baseName}.distance_filtered.csv", 'w', newline = '') as csv_out, \
        open("${fasta.baseName}.distance_filtered.fa", 'w') as fasta_out:

        reader = csv.DictReader(csv_in, delimiter=",", quotechar='\"', dialect = "unix")
        writer = csv.DictWriter(csv_out, fieldnames = reader.fieldnames, delimiter=",", quotechar='\"', quoting=csv.QUOTE_MINIMAL, dialect = "unix")
        writer.writeheader()

        for row in reader:
            if row["why_excluded"]:
                writer.writerow(row)
                continue
            id = row["sequence_name"]
            if id in reject:
                row["why_excluded"] = "distance to WH04 more than 4.0 epi-week std devs"
                writer.writerow(row)
                continue
            if id in alignment:
                writer.writerow(row)
                seq = str(alignment[id].seq)
                fasta_out.write(">" + id + "\\n")
                fasta_out.write(seq + "\\n")
    """
}


workflow postprocess_gisaid {
    take:
        gisaid_fasta
        gisaid_metadata
    main:
        gisaid_counts_by_country(gisaid_metadata)
        gisaid_distance_QC(gisaid_fasta, gisaid_metadata)
        gisaid_filter_on_distance_to_WH04(gisaid_fasta, gisaid_metadata,gisaid_distance_QC.out)
    emit:
        fasta = gisaid_filter_on_distance_to_WH04.out.fasta
        metadata = gisaid_filter_on_distance_to_WH04.out.metadata
}


workflow {
    gisaid_fasta = file(params.gisaid_fasta)
    gisaid_metadata = file(params.gisaid_metadata)

    postprocess_gisaid(gisaid_fasta, gisaid_metadata)
}
