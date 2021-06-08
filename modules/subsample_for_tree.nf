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
    path "${fasta.baseName}.site_filtered.fa"

    script:
    if ( params.downsample )
        """
        $project_dir/../bin/filter_by_ambiguous_sites.py \
                --in-alignment ${fasta} \
                --out-alignment "${fasta.baseName}.site_filtered.fa" \
                --sites ${ambiguous_sites}
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
        downsample(filter_on_ambiguous_sites.out, filter_on_sample_date.out)
        announce_summary(fasta, hash_non_unique_seqs.out.fasta, filter_on_ambiguous_sites.out, filter_on_sample_date.out, downsample.out.fasta)
    emit:
        fasta = downsample.out.fasta // subset of unique
        metadata = downsample.out.metadata
        hashmap = hash_non_unique_seqs.out.hashmap
}



workflow {
    fasta = file(params.fasta)
    metadata = file(params.metadata)

    subsample_for_tree(fasta,metadata)
}
