# Selection pipeline

From CDS orthologs to selection tests

## Intro

Each step of the pipeline is contained within one script, sometimes python sometimes snakemake. For python, some paths may be hard coded. For snakemake, a config file () provides all the paths and will need to be updated for your file system. There is also a snakemake profile () that allocates resources for your given cluster.

## 0. Setup

### Expected input files

Unaligned orthologous coding sequences (codons) from each genome. While your ortho predictions are probably done with peptide sequences, you will need to retrieve the corresponding CDS sequences to use this pipeline. I can provide advice for this part as well but it really depends on how the data is formatted at the beginning and whether or not you've done ortholog calling yet.

### Pipeline conda environment

I have provided a conda environment file `environment.yml` that contains all the software required for the pipeline, notably [Guidance](http://guidance.tau.ac.il/overview), [IQ-TREE](http://www.iqtree.org/), [HypPhy](http://www.hyphy.org/), and [snakemake](https://snakemake.readthedocs.io/en/stable/) and their dependencies. To install the environment with `conda`, run the following. 

```bash
conda env create -f environment.yml
```

This will create the environment as specifed in the file. The environment name is `cds-aln`, so to activate the environment, type:

```bash
conda activate cds-aln
```

If you don't have conda installed, see https://docs.conda.io/en/latest/miniconda.html

Since this is a rather large environment, and `conda` can be slow, I also recommend installing and using `mamba` instead of `conda`. This is a faster version of `conda`: https://mamba.readthedocs.io/en/latest/installation.html

If you try this, just replace the word `conda` in the above commands with `mamba`.

With the environment activated, try running the following commands: 

```bash
guidance
```

```bash
iqtree
```

```bash
hyphy -h
```

```bash
snakemake -h
```

In each case, if the program is properly installed, a help menu should come up. If you get an error instead of the help menu for any of these commands, something didn't install properly.

### The snakemake config file

The snakemake config file is a file that the snakemake scripts will read and should have the paths to your data in it. It is located at `scripts/snakemake-config.yaml` and currently has the paths to the files/directories on my file system.

You will have to edit the paths for your file system. 

In most cases, these are directories that **will be created by snakemake** during a given step, so they need not exist yet. The only one that must exist in the beginning is the `cds_directory` that contains FASTA files with CDS sequences (codons) for each orthogroup to be processed.

The `aln_filter_file` won't exist until after step 2 and will need to be filled in after that.

### The snakemake profile

Snakemake can submit jobs to a compute cluster, but it needs to know the rules and resources of the cluster. That is what `scripts/slurm_profile/config.yaml` is for. Set the `default-resources` here to be specific for your cluster and needs and add any other options you need. Things like `{rule}` and `{wildcards}` are inherited from the snakemake script itself.

The `jobs` parameter specifies how many individual jobs are submitted at once.

This is obviously for a SLURM cluster. If you have something else, let me know.

## 1. Align orthologous CDS sequences with Guidance

With the CDS sequences for each orthogroup in a directory, set that directory as `cds_directory` in the `scripts/snakemake-config.yaml` file. Set a desired directory as the output directory for the alignments as `aln_directory`. This will be created by snakemake if it doesn't exist.

At this point, I would enter the scripts directory with `cd scripts` and execute a dry run of `01_guidance.smk`:

```bash
snakemake -p -s 01_guidance.smk --configfile snakemake-config.yaml --profile slurm_profile/ --keep-going --dryrun
```

<details><summary>Click here to see a breakdown the above command</summary>
<p>

| Command line parameter  | Description |
|-------------------------|-------------|
| snakemake                          | A workflow management program we will use to submit and run jobs |
| -p                                 | Print out the shell commands that will be executed |
| -s 01_guidance.smk                 | The snakefile with rules (jobs) to be executed |
| --configfile snakemake-config.yaml | The config file read by the snakefile with information specific to our analysis |
| --profile slurm_profile/           | The directory containing the cluster profile (`config.yaml`) used by snakemake to allocate resources and distribute jobs to a compute cluster  |
| --keep-going                       | This tells snakemake to continue with other jobs if one or more fail. With such a large number of jobs some are likely to fail for some reason or another and this ensures that it doesn't hold up the rest of the jobs. The failed jobs should be looked at after others have finished |
| --dryrun                           | This tells snakemake to analyze the input and determine which jobs need to be run, but not to run them |

</p>
</details>
</br>

If the dry run completes successfully, you will quickly see every command that is to be run printed to the screen (thanks to `-p`), a list of all output files that will be generated in green, and then something like this:

```
Job stats:
job             count    min threads    max threads
------------  -------  -------------  -------------
all                 1              1              1
run_guidance     3052              1              1
total            3053              1              1


This was a dry-run (flag -n). The order of jobs does not reflect the order of execution.
```

If everything looks in order at this point, re-run the above command but remove the `--dryrun` option. This will submit each job individually.

Note that I highly recommend that when you submit the snakemake script after the dry run, you do this from a connection that is unlikely to be interrupted, either by submitting it as its own job on your cluster, or by activating a [screen](https://en.wikipedia.org/wiki/GNU_Screen).

Afterwards, address any failed jobs. These can be identified by checking the log files in the `aln_directory` `logs/01-Guidance/` subfolder (e.g. `grep "error\|Error\|ERROR" *.log`).

## 2. Filter alignments

For more confident analyses of selection, alignment filtering must be done to remove gappy and poorly aligned sites. To do this, I have a python script called `02_aln_filter.py`. 

Note that Guidance itself assigns a confidence score to each column in the alignment and filters out low quality sites. This is the main reason I use Guidance, and you should read the docs to see more info about these settings (which I just left at default). However, we also want to filter for things like overall gappiness of an alignment, and putatively mis-aligned sections of the alignment. `02_aln_filter.py` does this by:

1. Removing sequences in an alignment that have a percentage of gaps over some thresholed (`-s`, default 20%)
2. Removing alignment columns based on a sliding window of three codons: if `-c` percent of sequences for a given window have 2 codons with at least 2 gaps, then that 3 codon window is removed from all sequences (default 50%).
3. Removing alignments with fewer than 4 species and removing alignments shorter than 33 codons long.

Run this script as follows:

```bash
python 02_aln_filter.py -i [directory with Guidance alignments] -o [desired output directory for filtered alignments and filter summary files]
```

<details><summary>Click here to see a breakdown the above command</summary>
<p>

| Command line parameter  | Description |
|-------------------------|-------------|
| python 02_aln_filter.py                                                        | Run the script with Python |
| -i [directory with Guidance alignments]                                        | The directory containing the Guidance alignments produced from step 1. This will be the same as the `aln_directory` specified in the `snakemake-config.yaml` file. |
| -o [desired output directory for filtered alignments and filter summary files] | The output directory. Sequences will be placed in a sub-directory called `cds/` |

</p>
</details>
</br>

This will take several hours to run.

## 3. Making trees (optional if species tree already obtained)

At this point, you'll need to consider whether phylogenetic discordance is a major concern in your dataset. If it is, you will want to run selection tests on individual gene trees. However, as [we found with the first run](https://gwct.github.io/fish-genomics/fish_genomics.html#species-tree-and-gene-concordance-factors), with such deep divergence, discordance may not be an issue, and if you already have a species tree it may be much faster just to use that. 

However, if you do this, to assuage possible reviewer concerns, I would recommend demonstrating that discordance is not prevalent in the data. The easiest way to do this would be with [site concordance factors in IQ-Tree](http://www.iqtree.org/doc/Concordance-Factor). This would provide a concrete quantification of discordance and obviate the need for making gene trees.

If phylogenetic discordance is a concern, or if you do not have a species tree, you will have to generate gene trees and a species tree.

### Gene trees

I have provided `03_iqtree.smk` to make gene trees from the filtered alignments. First, update `snakemake-config.yaml` to make sure it has the correct `aln_filter_directory` and desired `tree_directory`, which will be created when this step is run. Then, executed a dry run as follows:

```bash
snakemake -p -s 03_iqtree.smk --configfile snakemake-config.yaml --profile slurm_profile/ --keep-going --dryrun
```

<details><summary>Click here to see a breakdown the above command</summary>
<p>

| Command line parameter  | Description |
|-------------------------|-------------|
| snakemake                          | A workflow management program we will use to submit and run jobs |
| -p                                 | Print out the shell commands that will be executed |
| -s 03_iqtree.smk                   | The snakefile with rules (jobs) to be executed |
| --configfile snakemake-config.yaml | The config file read by the snakefile with information specific to our analysis |
| --profile slurm_profile/           | The directory containing the cluster profile (`config.yaml`) used by snakemake to allocate resources and distribute jobs to a compute cluster  |
| --keep-going                       | This tells snakemake to continue with other jobs if one or more fail. With such a large number of jobs some are likely to fail for some reason or another and this ensures that it doesn't hold up the rest of the jobs. The failed jobs should be looked at after others have finished |
| --dryrun                           | This tells snakemake to analyze the input and determine which jobs need to be run, but not to run them |

</p>
</details>
</br>

Again, if the dry run looks ok, go ahead and remove that option and execute the full workflow to create gene trees for all alignments. Afterwards, address any failed jobs. These can be identified by checking the log files in the `tree_directory` `logs/cds-iqtree/` subfolder (e.g. `grep "error\|Error\|ERROR" *.log`).

### Species tree

A species tree could be created by concatenating the input alignments with IQ-Tree or by summarizing the gene trees with ASTRAL. Let me know if you need further info about these steps.

## 4. Selection tests

Now, we can analyze the aligned coding sequences to determine which lineages have experienced accelerated evolution in each gene. To do this, I've implemented a snakemake rule to run [HyPhy's aBSREL test](http://www.hyphy.org/methods/selection-methods/), which provides a p-value for each branch in a tree for each gene. Other tests are available (see linked page) and can easily be implemented in a similar fashion.

This rule is located in `04_hyphy.smk`. Be sure to update `snakemake-config.yaml` to include the correct path to the `aln_filter_directory` and the `species_tree_file`. If running with gene trees instead, we will have to modify the snakefile to use `tree_directory`. Contact me if this is the case. Otherwise, execute a dry run of the rule:

```bash
snakemake -p -s 04_hyphy.smk --configfile snakemake-config.yaml --profile slurm_profile/ --keep-going --dryrun
```

<details><summary>Click here to see a breakdown the above command</summary>
<p>

| Command line parameter  | Description |
|-------------------------|-------------|
| snakemake                          | A workflow management program we will use to submit and run jobs |
| -p                                 | Print out the shell commands that will be executed |
| -s 04_hyphy.smk                    | The snakefile with rules (jobs) to be executed |
| --configfile snakemake-config.yaml | The config file read by the snakefile with information specific to our analysis |
| --profile slurm_profile/           | The directory containing the cluster profile (`config.yaml`) used by snakemake to allocate resources and distribute jobs to a compute cluster  |
| --keep-going                       | This tells snakemake to continue with other jobs if one or more fail. With such a large number of jobs some are likely to fail for some reason or another and this ensures that it doesn't hold up the rest of the jobs. The failed jobs should be looked at after others have finished |
| --dryrun                           | This tells snakemake to analyze the input and determine which jobs need to be run, but not to run them |

</p>
</details>
</br>

If everything looks ok, remove the dry run option and execute! Afterwards, address any failed jobs. These can be identified by checking the log files in the `hyphy_directory` `logs/absrel/` subfolder (e.g. `grep "error\|Error\|ERROR" *.log`).

## 5. Post-analysis

When we get to this point, if you need help analyzing the results let me know. I have a [separate repository](https://github.com/gwct/hyphy-interface) for organizing the results of HyPhy runs that can be used.