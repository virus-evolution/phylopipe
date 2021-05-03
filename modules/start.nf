#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
publish_dev = file(params.publish_dev)

process get_git_hash {
    /**
    * Gets git commit
    */
    publishDir "${publish_dev}", mode: 'copy', overwrite: true

    input:
    path commit_file

    output:
    path "${commit_file}"

    script:
    """
    echo "\n Git hash \t = \t \$( git rev-parse HEAD) \n\n" >> ${commit_file}
    """
}

workflow start {
    params_file = file("${workDir}/input_params.txt")
    params_file << "\n#######################################################################################\n\n"

    printMapClosure = { key, value ->
            params_file << "$key = $value\n"
    }
    params.each(printMapClosure)
    get_git_hash(params_file)
}
