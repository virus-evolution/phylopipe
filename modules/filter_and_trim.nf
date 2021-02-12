#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

project_dir = projectDir
publish_dev = file(params.publish_dev)


process filter_low_coverage_sequences {
    /**
    * Keeps only sequences with completeness greater than min_covg threshold
    * @input alignment, metadata
    * @output alignment_updated, metadata_updated
    * @params min_covg
    */

    publishDir "${publish_dev}/${params.category}", pattern: "*.csv", mode: 'copy', saveAs: {"${params.category}_master.csv"}

    input:
    path alignment
    path metadata

    output:
    path "${alignment.baseName}.low_covg_filtered.fasta", emit: fasta_updated
    path "${metadata.baseName}.low_covg_filtered.csv", emit: metadata_updated

    script:
    if (!params.min_covg)
        """
        mv "${alignment}" "${alignment.baseName}.low_covg_filtered.fasta"
        mv "${metadata}" "${metadata.baseName}.low_covg_filtered.csv"
        """
    else
        """
        #!/usr/bin/env python3
        from Bio import SeqIO
        import csv

        alignment = SeqIO.index("${alignment}", "fasta")

        with open("${metadata}", 'r', newline = '') as csv_in, \
             open("${metadata.baseName}.low_covg_filtered.csv", 'w', newline = '') as csv_out, \
             open("${alignment.baseName}.low_covg_filtered.fasta", 'w') as fasta_out:

            reader = csv.DictReader(csv_in, delimiter=",", quotechar='\"', dialect = "unix")
            writer = csv.DictWriter(csv_out, fieldnames = reader.fieldnames, delimiter=",", quotechar='\"', quoting=csv.QUOTE_MINIMAL, dialect = "unix")
            writer.writeheader()

            for row in reader:
                if row["why_excluded"]:
                    writer.writerow(row)
                    continue
                id = row["sequence_name"]
                if id in alignment:
                    seq = str(alignment[id].seq)
                    mapped_completeness = float(len(seq.replace("N", "")) / len(seq))
                    if mapped_completeness >= float(${params.min_covg} / 100):
                        writer.writerow(row)
                        fasta_out.write(">" + id + "\\n")
                        fasta_out.write(seq + "\\n")
                    else:
                        row["why_excluded"] = "low mapped_completeness"
                        writer.writerow(row)
        """
}


process trim_alignment {
    /**
    * Trims start and end of alignment
    * @input alignment
    * @output alignment_updated
    * @params trim_start, trim_end
    */

    input:
    path alignment

    output:
    path "${alignment.baseName}.trimmed.fa"

    script:
    if (params.trim_start && params.trim_end)
        """
        #!/usr/bin/env python3
        from Bio import SeqIO

        strt = int(${params.trim_start})
        stp = int(${params.trim_end})

        with open("${alignment}", "r") as fasta_in, \
             open("${alignment.baseName}.trimmed.fa", "w") as fasta_out:

            for record in SeqIO.parse(fasta_in, "fasta"):
                seq = str(record.seq).upper()
                new_seq = ("N" * strt) + seq[strt:stp] + ("N" * (len(seq) - stp))
                fasta_out.write(">" + record.id + "\\n")
                fasta_out.write(new_seq + "\\n")
        """
    else
        """
        mv "${alignment.baseName}" "${alignment.baseName}.trimmed.fa"
        """
}


workflow filter_and_trim {
    take:
        in_fasta
        in_metadata
    main:
        filter_low_coverage_sequences(in_fasta, in_metadata)
        trim_alignment(filter_low_coverage_sequences.out.fasta_updated)
    emit:
        fasta = trim_alignment.out
        metadata = filter_low_coverage_sequences.out.metadata_updated
}

workflow {
    in_fasta = file(params.gisaid_fasta)
    in_metadata = file(params.gisaid_metadata)

    filter_and_trim(in_fasta,in_metadata)
}
