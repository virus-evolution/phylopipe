#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

project_dir = projectDir


process gisaid_process_json {
    /**
    * Downloads
    * @input json
    * @output fasta, metadata
    * @params omissions
    */

    input:
    path json

    output:
    path "gisaid.fasta", emit: fasta
    path "gisaid.csv", emit: metadata

    script:
    """
    datafunk process_gisaid_data \
        --input-json ${json} \
        --input-metadata False \
        --exclude-file ${omissions} \
        --output-fasta "gisaid.fasta" \
        --output-metadata "gisaid.csv" \
        --exclude-undated
    """
}


process gisaid_add_columns_to_metadata {
    input:
    path gisaid_fasta
    path gisaid_metadata

    output:
    path "${gisaid_metadata.baseName}.add_metadata.csv"

    script:
    """
    #!/usr/bin/env python3
    from Bio import SeqIO
    import csv

    alignment = SeqIO.index("${fasta}", "fasta")

    with open("${gisaid_metadata}", 'r', newline = '') as csv_in, \
        open("${gisaid_metadata.baseName}.add_metadata.csv", 'w', newline = '') as csv_out:

        reader = csv.DictReader(csv_in, delimiter=",", quotechar='\"', dialect = "unix")
        writer = csv.DictWriter(csv_out, fieldnames = reader.fieldnames + ['sequence_name', 'why_excluded'], delimiter=",", quotechar='\"', quoting=csv.QUOTE_MINIMAL, dialect = "unix")
        writer.writeheader()

        for row in reader:
            edin_header = row["edin_header"]
            new_header = edin_header.split("|")[0]
            row['sequence_name'] = new_header
            if edin_header not in alignment:
                row['why_excluded'] = "filtered during loading from JSON"
            else:
                row['why_excluded'] = ""
            writer.writerow(row)
    """
}


omissions = file(params.omissions)

workflow preprocess_gisaid {
    take:
        gisaid_json
    main:
        gisaid_process_json(gisaid_json)
        gisaid_add_columns_to_metadata(gisaid_process_json.out.metadata)
    emit:
        fasta = gisaid_process_json.out.fasta
        metadata = gisaid_add_columns_to_metadata.out
}


workflow {
    gisaid_json = file(params.gisaid_json)

    preprocess_gisaid(gisaid_json)
}
