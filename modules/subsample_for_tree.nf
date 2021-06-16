#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

project_dir = projectDir
publish_dev = file(params.publish_dev)


process hash_non_unique_seqs {
    /**
    * Subsets a unique set of sequences
    * @input fasta, metadata
    */
    memory { fasta.size() * 2.B  + 2.GB * task.attempt }
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
        $project_dir/../bin/date_filter.py \
                    --in-metadata ${metadata} \
                    --out-metadata "${metadata.baseName}.date_filtered.csv" \
                    --date ${params.date} \
                    --time-window ${params.time_window} \
                    --filter-column "date_filter"
        """
    else
        """
        echo "date ${params.date}"
        echo "time window ${params.time_window}"
        mv "${metadata}" "${metadata.baseName}.date_filtered.csv"
        """
}

process filter_on_ambiguous_sites {
    /**
    * Filter to keep only sequences without ambiguous bases at those sites
    * @input fasta
    * @output fasta
    * @params sites
    */

    input:
    path fasta

    output:
    path "${fasta.baseName}.site_filtered.fa", emit: fasta
    path "ids_with_ambiguous_sites.log", emit: ambiguous_site_ids


    script:
    if ( params.downsample )
        """
        $project_dir/../bin/filter_by_ambiguous_sites.py \
                --in-alignment ${fasta} \
                --out-alignment "${fasta.baseName}.site_filtered.fa" \
                --sites ${ambiguous_sites} > "ids_with_ambiguous_sites.log"
        """
    else
        """
        mv "${fasta}" "${fasta.baseName}.site_filtered.fa"
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

process annotate_metadata {
    /**
    * Adds note column with info about sequences hashed, filtered and downsampled
    * @input metadata
    * @output metadata
    */

    input:
    path metadata
    path hashmap
    path ambiguous_site_log

    output:
    path "${metadata.baseName}.annotated.csv"

    script:
    """
    #!/usr/bin/env python3
    import csv

    tips = set()
    hashmap = {}
    with open("${hashmap}", 'r', newline = '') as hashmap_in:
        for line in hashmap_in:
            tip, redundant = line.rstrip().split(",")
            tips.add(tip)
            for id in redundant.split("|"):
                hashmap[id] = tip
    print("hashmap contains %d tips and %d redundants" %(len(tips), len(hashmap)))

    ambig = set()
    with open("${ambiguous_site_log}", 'r', newline = '') as ambiguous_site_log_in:
        for line in ambiguous_site_log_in:
            ambig.add(line.rstrip())
    print("%d had ambiguous bases" %len(ambig))

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
            if row["sequence_name"] in hashmap:
                note.append("hashed to tip %s" %hashmap[row["sequence_name"]])
            if row["sequence_name"] in ambig:
                note.append("filtered due to ambiguous base")
            if row["date_filter"]:
                if row["note"] and "downsample" not in row["note"]:
                    note.append("date filtered")
            if row["note"]:
                note.append(row["note"])
            row["note"] = "|".join(note)
            writer.writerow(row)
    """
}

process announce_summary {
    /**
    * Summarizes subsampling into JSON
    * @input fasta
    */

    input:
    path original
    path unique
    path filtered_on_ambiguous_sites
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
            echo "> Number of unique sequences : \$(cat ${unique} | grep '>' | wc -l)\\n" >> announce.json
            echo "> Number of sequences with non-ambiguous bases at sites of interest : \$(cat ${filtered_on_ambiguous_sites} | grep '>' | wc -l)\\n" >> announce.json
            echo "> Number of (non-unique) sequences with sample_date older than ${params.time_window} days: \$(cat ${filtered_on_sample_date} | grep 'sample_date older than' | wc -l)\\n" >> announce.json
            echo "> Number of sequences after downsampling: \$(cat ${downsampled} | grep '>' | wc -l)\\n" >> announce.json
            echo '"}' >> announce.json

            echo 'webhook ${params.webhook}'

            curl -X POST -H "Content-type: application/json" -d @announce.json ${params.webhook}
            """
        else
            """
            echo '{"text":"' > announce.json
            echo "*${params.whoami}: Subsampling ${params.date} for tree*\\n" >> announce.json
            echo "> Number of sequences in COG and GISAID input files : \$(cat ${original} | grep '>' | wc -l)\\n" >> announce.json
            echo "> Number of unique sequences : \$(cat ${unique} | grep '>' | wc -l)\\n" >> announce.json
            echo "> Number of sequences with non-ambiguous bases at sites of interest : \$(cat ${filtered_on_ambiguous_sites} | grep '>' | wc -l)\\n" >> announce.json
            echo "> Number of (non-unique) sequences with sample_date older than ${params.time_window} days: \$(cat ${filtered_on_sample_date} | grep 'sample_date older than' | wc -l)\\n" >> announce.json
            echo "> Number of sequences after downsampling: \$(cat ${downsampled} | grep '>' | wc -l)\\n" >> announce.json
            echo '"}' >> announce.json
            """
}

lineage_splits = file(params.lineage_splits)
mask_file = file(params.mask)
ambiguous_sites = file(params.ambiguous_sites)


workflow subsample_for_tree {
    take:
        fasta
        metadata
    main:
        hash_non_unique_seqs(fasta, metadata)
        filter_on_ambiguous_sites(hash_non_unique_seqs.out.fasta)
        filter_on_sample_date(metadata)
        downsample(filter_on_ambiguous_sites.out.fasta, filter_on_sample_date.out)
        annotate_metadata(downsample.out.metadata, hash_non_unique_seqs.out.hashmap, filter_on_ambiguous_sites.out.ambiguous_site_ids)
        announce_summary(fasta, hash_non_unique_seqs.out.fasta, filter_on_ambiguous_sites.out.fasta, filter_on_sample_date.out, downsample.out.fasta)
    emit:
        fasta = downsample.out.fasta // subset of unique
        metadata = annotate_metadata.out
        hashmap = hash_non_unique_seqs.out.hashmap
}



workflow {
    fasta = file(params.fasta)
    metadata = file(params.metadata)

    subsample_for_tree(fasta,metadata)
}
