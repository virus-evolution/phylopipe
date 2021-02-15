#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

project_dir = projectDir
publish_dev = file(params.publish_dev)


process build_start_tree {
    /**
    * Makes an initial tree with fasttree
    * @input alignment
    */

    input:
    path alignment

    output:
    path "${alignment}.1.multi.fasttree"

    script:
    """
    fasttree -nt -gamma -nosupport \
      -sprlength 500 \
      -nni 0 \
      -spr 5 \
      -refresh 0.8 \
      -topm 1.5 \
      -close 0.75 \
      -noml ${alignment} > "${alignment}.multi.fasttree"

    fasttree -nt -gamma \
      -sprlength 200 \
      -spr 5 \
      -intree "${alignment}.multi.fasttree" ${alignment} > "ft_SH.unrooted.tree"
    """
}


process bootstrap_trees {
    /**
    * Makes bootstrap trees
    * @input val, alignment
    */

    input:
    tuple path(alignment), val(x)

    output:
    path "bootstrap_unrooted_${x}.tree"

    script:
    """
    goalign build seqboot -i ${alignment} -t 1 -n 1 -S -o bootstrap_${x}
    fasttree -nosupport -nt -fastest "bootstrap_${x}.fa" > "bootstrap_unrooted_${x}.tree"
    """
}


process get_support {
    /**
    * Downloads
    * @input first_tree, replicate_trees
    */

    input:
    path first_tree
    path replicate_trees

    output:
    path "ft_TBE.unrooted.tree", emit: tbe_tree
    path "ft_FBP.unrooted.tree", emit: fbp_tree

    script:
    """
    gotree compute support tbe -i ${first_tree} -b ${replicate_trees} -t $threads -o ft_TBE.unrooted.tree
    gotree compute support fbp -i ${first_tree} -b ${replicate_trees} -t $threads -o ft_FBP.unrooted.tree
    """
}


process reroot_tree {
    /**
    *
    * @input tree
    */

    input:
    path tree

    output:
    path "${tree.baseName}.rooted.tree"

    script:
    """
    nw_reroot ${tree} 'Wuhan/WH04/2020' > "${tree.baseName}.rooted.tree"
    """
}


process treeshrink {
    /**
    *
    * @input tree
    */

    input:
    path tree

    output:
    path "treeshrink/*.tree"

    script:
    """
    run_treeshrink.py -t ${tree} -q 0.05 -c -o treeshrink
    """
}

process reroot_treeshrink {
    /**
    * Downloads
    * @input tree
    */

    input:
    path tree

    output:
    path ${tree.baseName}.tree
    path alignments.log

    script:
    """
    nw_reroot ${tree} 'Wuhan/WH04/2020' > "${tree.baseName}.tree"
    sed -i.bak "s/'//g" "${tree.baseName}.tree"
    echo "Tree stats for ${tree.baseName}.tree" > alignments.log
    nw_stats "${tree.baseName}.tree" >> alignments.log
    """
}



workflow build_tree {
    take:
        alignment
    main:
        build_start_tree(alignment)
        bootstraps = Channel.from(1..100)
                            .combine(alignment)
        bootstrap_trees(bootstraps)
        bootstraps.out.collectFile(name: 'ft_replicates_multi.tree', newLine: true)
                      .set{ replicates }
        get_support(build_start_tree.out, replicates)
        build_start_tree.out.concat(get_support.out.tbe_tree, get_support.out.fbp_tree).set{ trees }
        reroot_tree(trees)
        treeshrink(reroot_tree.out)
        reroot_treeshrink(treeshrink.out)
    emit:
        reroot_treeshrink.out

}


workflow {
    alignment = file(params.alignment)

    build_tree(alignment)
}