#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

project_dir = projectDir
publish_dev = file(params.publish_dev)


process filter_on_sample_date_for_recent {
    /**
    * If a time window (in days) is provided, excludes samples from
    * METADATA files which do not fall within X days of date
    * @input metadata
    * @output metadata_updated
    * @params time_window, date
    */

    input:
    path metadata

    output:
    path "${metadata.baseName}.recent.csv"

    script:
    if ( params.time_window && params.date )
        """
        $project_dir/../bin/date_filter.py \
                    --in-metadata ${metadata} \
                    --out-metadata "${metadata.baseName}.recent.csv" \
                    --date ${params.date} \
                    --time-window ${params.time_window} \
                    --filter-column "why_excluded" \
                    --restrict
        """
    else
        """
        echo "date ${params.date}"
        echo "time window ${params.time_window}"
        mv "${metadata}" "${metadata.baseName}.recent.csv"
        """
}

process annotate_metadata {
    /**
    * Adds to note column with info about sequences excluded by date
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

    print("Fasta of recent had %d sequences" %len(seqs))

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
            if row["note"]:
                note.extend(row["note"].split("|"))
            statement = "sample not recent"
            if row["sequence_name"] not in seqs and statement not in note:
                note.append("sample not recent")
            row["note"] = "|".join(note)
            writer.writerow(row)
    """
}

process fetch_recent {
    /**
    * Fetches fasta of recent sequences
    * @input fasta, metadata
    * @output fasta
    */

    input:
    path fasta
    path metadata

    output:
    path "${fasta.baseName}.recent.fa"

    script:
    """
    fastafunk fetch \
        --in-fasta ${fasta} \
        --in-metadata ${metadata} \
        --out-fasta "${fasta.baseName}.recent.fa" \
        --index-column "sequence_name" \
        --low-memory
    """
}


process fetch_designations {
    /**
    * Extracts fasta corresponding to lineage designations fasta
    * @input fasta
    */
    label 'retry_increasing_mem'

    input:
    path fasta
    path metadata

    output:
    path "${fasta.baseName}.designations.fasta"

    script:
    """
    fastafunk extract \
        --in-fasta ${fasta} \
        --in-metadata ${metadata} \
        --out-fasta "${fasta.baseName}.designations.fasta" \
        --reject-fasta "tmp.fasta" \
        --low-memory
    """
}


process fetch_outgroups {
    /**
    * Extracts fasta corresponding to lineage outgroups fasta
    * @input fasta
    */
    label 'retry_increasing_mem'

    input:
    path fasta

    output:
    path "${fasta.baseName}.outgroups.fasta"

    script:
    """
    tail -n+2 ${lineage_splits} | cut -f2 -d"," >> outgroups.csv
    fastafunk extract \
        --in-fasta ${fasta} \
        --in-metadata "outgroups.csv" \
        --out-fasta "${fasta.baseName}.outgroups.fasta" \
        --reject-fasta "tmp.fasta" \
        --low-memory
    """
}


lineage_splits = file(params.lineage_splits)


workflow extract_protected_sequences {
    take:
        fasta
        metadata
    main:
        if ( params.time_window && params.date ) {
            filter_on_sample_date_for_recent(metadata)
            fetch_recent(fasta, filter_on_sample_date_for_recent.out)
            recent_ch = fetch_recent.out
        } else {
            recent_ch = Channel.empty()
        }

        if ( params.lineage_designations ) {
            lineage_designations = file(params.lineage_designations)
            fetch_designations(fasta, lineage_designations)
            designations_ch = fetch_designations.out
        } else {
            designations_ch = Channel.empty()
        }

        fetch_outgroups(fasta)
        fetch_outgroups.out.concat( recent_ch, designations_ch ).set { protected_ch }
        protected_ch.collectFile(name: "force.fa").set { output_ch }
        annotate_metadata(metadata,output_ch)
    emit:
        fasta = output_ch
        metadata = annotate_metadata.out
}



workflow {
    fasta = file(params.fasta)
    metadata = file(params.metadata)

    extract_protected_sequences(fasta, metadata)
}
