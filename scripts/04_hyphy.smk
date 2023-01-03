#############################################################################
# Snakemake rule to run HyPhy's aBSREL on an input directory with CDS alignments
# Gregg Thomas, January 2023
#############################################################################

# snakemake -p -s 01_guidance.smk --configfile snakemake-config.yaml --profile slurm_profile/ --keep-going --dryrun

#############################################################################

import os

#############################################################################

ALNDIR = config["aln_filter_directory"]
TREEDIR = config["tree_directory"]
SPECTREE = config["species_tree_file"]
HYPHYDIR = config["hyphy_directory"]

loci = [ cds_file.split("-cds.guidance.filter.fa")[0] for cds_file in os.listdir(CDSDIR) ];
print("# running selection tests on ", len(loci), " loci");

#############################################################################
# Final rule - rule that depends on final expected output file and initiates all
# the other rules

localrules: all

rule all:
    input:
        expand(os.path.join(HYPHYDIR, "absrel", "{locus}-cds.json"), locus=loci)
        # Expected output from rule absrel

#############################################################################
# Pipeline rules

rule absrel:
    input:
        aln = os.path.join(ALNDIR, "02-Filter-spec7-seq50-site50", "cds-spec", "{locus}-cds.guidance.filter.fa"),
        tree = SPECTREE
    output:
        os.path.join(HYPHYDIR, "absrel", "{locus}-cds.json")
    log:
        os.path.join(HYPHYDIR, "logs", "absrel", "{locus}-cds.log")
    shell:
        """
        hyphy absrel --alignment {input.aln} --tree {input.tree} --output {output} &> {log}
        """
# Run each locus through absrel

#############################################################################