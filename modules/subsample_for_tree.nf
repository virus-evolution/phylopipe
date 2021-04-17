#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

project_dir = projectDir
publish_dev = file(params.publish_dev)


process mask_alignment {
    /**
    * Applies a mask to aligned FASTA
    * @input alignment
    * @output alignment_updated
    * @params mask_file
    */

    input:
    path alignment

    output:
    path "${alignment.baseName}.masked.fa"

    script:
    """
    $project_dir/../bin/add_mask.py \
      --in-alignment ${alignment} \
      --out-alignment "${alignment.baseName}.masked.fa" \
      --mask ${mask_file} \
    """
}

process filter_uk {
    /**
    * Filters UK sequences by metadata to exclude biosample duplicates
    * @input fasta, metadata
    */

    input:
    path fasta
    path metadata

    output:
    path "${fasta.baseName}.filtered.fasta", emit: fasta
    path "${metadata.baseName}.filtered.csv", emit: metadata

    script:
    """
    $project_dir/../bin/filter.py \
            --in-fasta ${fasta} \
            --in-metadata ${metadata} \
            --outgroups ${lineage_splits} \
            --out-fasta "${fasta.baseName}.filtered.fasta" \
            --out-metadata "${metadata.baseName}.filtered.csv" \
            --exclude_true duplicate
    """
}


process hash_non_unique_seqs {
    /**
    * Subsets a unique set of sequences
    * @input fasta, metadata
    */
    memory { 25.GB  + 2.GB * task.attempt }
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries = 2

    input:
    path fasta
    path metadata

    output:
    path "${fasta.baseName}.unique.fasta", emit: fasta
    path "${fasta.baseName}.hashmap.csv", emit: hashmap

    script:
    """
    $project_dir/../bin/hash_non_unique_seqs.py \
        --in-fasta ${fasta} \
        --in-metadata ${metadata} \
        --out-fasta "${fasta.baseName}.unique.fasta" \
        --out-metadata "${fasta.baseName}.hashmap.csv" \
        --outgroups ${lineage_splits}
    """
}

process filter_on_sample_date {
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
    path "${metadata.baseName}.date_filtered.csv"

    script:
    if ( params.time_window && params.date)
        """
        #!/usr/bin/env python3
        import datetime
        import csv

        window = datetime.timedelta(int("${params.time_window}"))
        todays_date = datetime.datetime.strptime("${params.date}", '%Y-%m-%d').date()
        print(window, todays_date)

        with open("${metadata}", 'r', newline = '') as csv_in, \
            open("${metadata.baseName}.date_filtered.csv", 'w', newline = '') as csv_out:

            reader = csv.DictReader(csv_in, delimiter=",", quotechar='\"', dialect = "unix")
            fieldnames = reader.fieldnames
            if "why_excluded" not in reader.fieldnames:
                fieldnames.append("why_excluded")
            fieldnames.append("date_filter")
            writer = csv.DictWriter(csv_out, fieldnames = fieldnames, delimiter=",", quotechar='\"', quoting=csv.QUOTE_MINIMAL, dialect = "unix")
            writer.writeheader()

            for row in reader:
                row["date_filter"] = ""

                if "why_excluded" not in row:
                    row["why_excluded"] = ""

                if row["why_excluded"] in reader.fieldnames and row["why_excluded"] not in [None,"None",""]:
                    writer.writerow(row)
                    continue

                try:
                    date = datetime.datetime.strptime(row["sample_date"], '%Y-%m-%d').date()
                except:
                    row["why_excluded"] = "no sample_date"
                    writer.writerow(row)
                    continue

                if (todays_date - window) > date:
                    row["date_filter"] = "sample_date older than %s days" %window
                    writer.writerow(row)
                    continue

                writer.writerow(row)
        """
    else
        """
        mv "${metadata}" "${metadata.baseName}.date_filtered.csv"
        """
}


process downsample {
    /**
    * Selects a subsample of sequences outside of the date region to include for background
    * @input metadata
    */

    input:
    path fasta
    path metadata

    output:
    path "${fasta.baseName}.downsampled.fa", emit: fasta
    path "${metadata.baseName}.downsampled.csv", emit: metadata

    script:
    if ( params.downsample )
        """

        $project_dir/../bin/downsample.py \
            --in-metadata ${metadata} \
            --in-fasta ${fasta} \
            --out-metadata "${metadata.baseName}.downsampled.csv" \
            --out-fasta "${fasta.baseName}.downsampled.fa" \
            --outgroups ${lineage_splits} \
            --diff ${params.downsample_diff} \
            --downsample_date_excluded
        """
    else
        """
        mv "${fasta}" "${fasta.baseName}.downsampled.fa"
        mv "${metadata}" "${metadata.baseName}.downsampled.csv"
        """
}


process announce_summary {
    /**
    * Summarizes subsampling into JSON
    * @input fasta
    */

    input:
    path original
    path deduplicated
    path unique
    path filtered_on_sample_date
    path downsampled

    output:
    path "announce.json"

    script:
        if (params.webhook)
            """
            echo '{"text":"' > announce.json
                echo "*${params.whoami}: Subsampling ${params.date} for tree*\\n" >> announce.json
                echo "> Number of sequences in COG and GISAID input files : \$(cat ${original} | grep '>' | wc -l)\\n" >> announce.json
                echo "> Number of sequences after filtering uk sequences (deduplication by biosample id): \$(cat ${deduplicated} | grep '>' | wc -l)\\n" >> announce.json
                echo "> Number of unique sequences : \$(cat ${unique} | grep '>' | wc -l)\\n" >> announce.json
                echo "> Number of (non-unique) sequences with sample_date older than ${params.time_window} days: \$(cat ${filtered_on_sample_date} | grep 'sample_date older than' | wc -l)\\n" >> announce.json
                echo "> Number of sequences after downsampling: \$(cat ${downsampled} | grep '>' | wc -l)\\n" >> announce.json
                echo '"}' >> announce.json

            echo 'webhook ${params.webhook}'

            curl -X POST -H "Content-type: application/json" -d @announce.json ${params.webhook}
            """
        else
            """
            touch "announce.json"
            """
}

lineage_splits = file(params.lineage_splits)
mask_file = file(params.mask)

workflow mask_and_filter {
    take:
        fasta
        metadata
    main:
        mask_alignment(fasta)
        filter_uk(mask_alignment.out, metadata)
    emit:
        masked_deduped_fasta = filter_uk.out.fasta
        metadata = filter_uk.out.metadata
}

workflow subsample_for_tree {
    take:
        fasta
        metadata
    main:
        mask_alignment(fasta)
        filter_uk(mask_alignment.out, metadata)
        hash_non_unique_seqs(filter_uk.out.fasta, filter_uk.out.metadata)
        filter_on_sample_date(filter_uk.out.metadata)
        downsample(hash_non_unique_seqs.out.fasta, filter_on_sample_date.out)
        announce_summary(fasta, filter_uk.out.fasta, hash_non_unique_seqs.out.fasta, filter_on_sample_date.out, downsample.out.fasta)
    emit:
        masked_deduped_fasta = filter_uk.out.fasta
        fasta = downsample.out.fasta // subset of unique
        metadata = downsample.out.metadata
        hashmap = hash_non_unique_seqs.out.hashmap
        unique = hash_non_unique_seqs.out.fasta
}



workflow {
    fasta = file(params.fasta)
    metadata = file(params.metadata)

    subsample_for_tree(fasta,metadata)
}
