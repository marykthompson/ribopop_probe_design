# Snakemake workflow: ribopop_probe_design

This workflow designs oligos that hybridize to rRNA targets. First, a consensus
sequence is built for each target. Next, potential probe target sites are screened against
the genome and the transcriptome for non-rRNA matches. Finally, probes are selected
to cover the targets, while adhering to the specified design parameters and avoiding
the predicted non-specific targeting.

## Authors

* Mary Kay Thompson (@marykthompson)
* Maria Kiourlappou (@mkiourlappou)

## Usage

#### Step 1: Install the workflow

To simply download the workflow, go to the 'Clone or download' tab and download
the ZIP file.

### Step 2: Prepare the design parameter files

Prepare an output directory. Prepare a config.yml file matching the format of example_outdir/config.yml. You need to prepare the following simple csv files, matching the format of those in example_outdir/design_params/

1. seqs_and_anns.csv: This file contains the urls or absolute paths to the gtf file and the genome, cdna, and ncrna fasta files for each organism in the target set.

2. targets.csv: This file contains the names of the rRNA targets, their corresponding organisms,
and the name of fasta file containing their sequences. Note that it is fine to provide all targets
in the same fasta file or to provide them in separate fasta files. If the fasta files contain sequences not matching the provided IDs, only those matching the provided IDs will be included. If left blank, the fasta file with default to the ncrna.fa file provided for that organism. Special regions in the target transcripts can also be indicated in this file (1-based, closed interval). The transcript-specific indices will be converted to consensus sequence indices during the alignment steps.
    excluded_regions: Regions that should not overlap with probes.
    target_subregions: Split the target up into multiple regions for probe design. This feature was used to specify that an equal number of probes should be targeted to the left and right side
  of the Drosophila 28S, which is cleaved into two fragments.

3. params.csv: This file specifies the design parameters, such as min and max Tm, number of probes
for each target, etc. The excluded_regions and target_subregions can instead be specified in this file
as excluded_regions_consensus and target_subregions_consensus, but then the indices must be relative to the consensus sequence.

### Step 3: (Recommended) Install Conda

If you use conda, you will not need to pre-install each required package.
We recommend the miniconda distribution. Please see here for instructions
on how to install miniconda for your operating system.
https://docs.conda.io/en/latest/miniconda.html

### Step 4: Set up the conda environment.

Create the probe design environment.

    conda env create -f ribopop_probe_design.yaml

### Step 5: Run the pipeline.

Activate the probe design environment.

    conda activate ribopop_probe_design

Run the pipeline.

    snakemake --directory <your output directory>

If you're not using conda, you will need to install snakemake as well as all the other
software listed in envs/ribopop_probe_design.yaml in advance of running the workflow using your preferred installation method.

### Additional options ###

It may be desirable to run probe design again on a different set of targets without
rebuilding the transcript index used for off-target screening. This can be accomplished
by running a different version of the Snakefile and providing a different parameter
file called indices.csv. In this case, you must also provide the path to the fasta
containing your target sequences in the targets.csv file. See example in example_outdir/design_params_prebuilt/. You can leave the target_homology field blank for the organism(s)
that you are using for building the target consensus, but you must provide it for any additional organisms that you are screening candidate probes against if they are not
included in building the target consensus.

Run the pipeline using a pre-built transcript index:

    snakemake -s use_prebuilt_indices.smk --directory <your output directory>

### Other notes ###

The pipeline assumes that fasta files used for building the transcript index end in
a newline character (\n) as the Ensembl ones used here do. If yours do not, you should
add one before running the pipeline.
