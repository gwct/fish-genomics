#############################################################################
# Snakemake rule to run Guidance on an input directory with CDS sequences
# Gregg Thomas, January 2023
#############################################################################

# snakemake -p -s 01_guidance.smk --configfile snakemake-config.yaml --profile slurm_profile/ --keep-going --dryrun

#############################################################################

import os

#############################################################################

ALNDIR = config["aln_directory"]
TREEDIR = config["tree_directory"]
FILTERFILE = config["aln_filter_file"]

filtered_loci = [ line.strip() for line in open(FILTERFILE, "r") if line.strip() != "" ];
print("# making gene trees for ", len(filtered_loci), " loci");

#############################################################################
# Final rule - rule that depends on final expected output file and initiates all
# the other rules

localrules: all

rule all:
    input:
        expand(os.path.join(TREEDIR, "cds-iqtree", "loci", "{filtered_locus}", "{filtered_locus}.treefile"), filtered_locus=filtered_loci)
        # Output from make_gene_trees

#############################################################################
# Pipeline rules

rule make_gene_trees:
    input: 
        os.path.join(ALNDIR, "02-Filter-spec7-seq50-site50", "cds", "{filtered_locus}-cds.guidance.filter.fa")
    output:
        os.path.join(TREEDIR, "cds-iqtree", "loci", "{filtered_locus}", "{filtered_locus}.treefile")
    params:
        prefix = os.path.join(TREEDIR, "cds-iqtree", "loci", "{filtered_locus}", "{filtered_locus}")
    log:
        os.path.join(TREEDIR,"logs", "cds-iqtree", "{filtered_locus}-cds-iqtree.log")
    resources:
        cpus = 4
    shell:
     """
        iqtree -s {input} --prefix {params.prefix} -B 1000 -T 4 &> {log}
     """

#############################################################################