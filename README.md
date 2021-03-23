## Phylopipe

Pipeline to create phylogenetic trees for UK and global SARS-CoV-2 sequences and metadata, and publish matched subsets of annotated trees, FASTA sequences and metadata for groups with different access to sensitive data. 

Builds trees weekly, with daily updates.

### Pipeline Overview

#### Preprocessing

1. Applies mask to problematic sites in the alignment `mask.txt`
2. Filters UK sequences, excluding those that have been labelled as duplicate `source_id`,  i.e. same patient
3. Hash non-unique sequences, storing a FASTA with a representative for each unique sequence, and a hashmap from the representative to the now excluded IDs with identical sequence. When an outgroups file is provided, protects the outgroups, and otherwise choses a representative from the most recent `epi-week` to prevent filtering by date downstream
4. Filter on sample date, keeping all sequences from the last 120 days, keeping specified outgroups, and downsampling the remainder by excluding sequences within 3 mutations from an included sequence

#### Build full tree by splitting and grafting

1. Split FASTA based on known distinct clades specified in `lineage_splits.csv`
2. Build tree for each sub-FASTA using `FastTreeMP` and reroot on the clade-specific outgrip
3. Graft together the subtrees to make a complete tree
4. Expand the hashmap, inserting polytomies for each non-unique sequence

#### Add downsample-excluded sequences where possible

1. Using `usher` and `faToVcf`, take the filtered aligned FASTA from preprocessing step 2 and construct a mutation annotated tree based on the grafted tree, adding the missing samples in the process where possible

#### Post-process tree

1. Sort and collapse short branches < 0.000005
2. Annotate tree tips with `country`, `lineage` and `uk_lineage`
3. Infer deltrans with ancestral reconstruction and annotate
4. Merge and create new `uk_lineages` and annotate
5. Infer `phylotypes` for UK lineages and annotate

#### Publish tree outputs

1. Publishes subsets of FASTA and metadata CSV with the NEWICK or NEXUS tree as specified in `publish_recipes.json`

#### Daily updates

1. Adds new sequences found in preprocessing step 2 to the `usher` mutation annotated tree daily, with the full tree pipeline run weekly.

### What is grapevine?

`grapevine` (https://github.com/COG-UK/grapevine) was the name of the original pipeline which preprocessed, aligned and variant called sequences, made phylogenetic trees and more. As the number of sequences has grown the tree building steps take increasingly long to complete. `Datapipe` (`https://github.com/COG-UK/grapevine_nextflow`) was created to provided daily alignment and metadata processing. This pipeline takes the output of datapipe, constructs trees, annotates and publishes them.
