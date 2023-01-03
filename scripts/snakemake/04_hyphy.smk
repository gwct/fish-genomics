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
SPECTREE = config["species_tree_file"]
FILTERFILE = config["aln_filter_file"]
HYPHYDIR = config["hyphy_directory"]

filtered_loci = [ line.strip() for line in open(FILTERFILE, "r") if line.strip() != "" ];
print("# running selection tests on ", len(filtered_loci), " loci");

#############################################################################
# Final rule - rule that depends on final expected output file and initiates all
# the other rules

localrules: all

rule all:
    input:
        expand(os.path.join(HYPHYDIR, "absrel", "{filtered_locus}-cds.json"), filtered_locus=filtered_loci)
        # Expected output from rule absrel

#############################################################################
# Pipeline rules

rule absrel:
    input:
        aln = os.path.join(ALNDIR, "02-Filter-spec7-seq50-site50", "cds-spec", "{filtered_locus}-cds.guidance.filter.fa"),
        tree = SPECTREE
    output:
        os.path.join(HYPHYDIR, "absrel", "{filtered_locus}-cds.json")
    log:
        os.path.join(HYPHYDIR, "logs", "absrel", "{filtered_locus}-cds.log")
    shell:
        """
        hyphy absrel --alignment {input.aln} --tree {input.tree} --output {output} &> {log}
        """
# Run each locus through absrel

#############################################################################