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
                    --filter-column "why_excluded"
        """
    else
        """
        echo "date ${params.date}"
        echo "time window ${params.time_window}"
        mv "${metadata}" "${metadata.baseName}.recent.csv"
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
            fetch_designations(fasta, metadata)
            designations_ch = fetch_designations.out
        } else {
            designations_ch = Channel.empty()
        }

        fetch_outgroups(fasta)
        fetch_outgroups.out.concat( recent_ch, designations_ch ).set { protected_ch }
    emit:
        protected_ch
}



workflow {
    fasta = file(params.fasta)
    metadata = file(params.metadata)

    extract_protected_sequences(fasta, metadata)
}
