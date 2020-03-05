# Snakemake workflow: ribopop_probe_design
This workflow designs oligos that hybridize to rRNA targets. First, a consensus
sequence is built for each target. Next, potential probes are screened against
the genome and the transcriptome for non-rRNA matches. Finally, probes are selected
to cover the targets, while adhering to the specified design parameters and avoiding
the predicted non-specific targeting.
## Authors

* Mary Thompson (@marykthompson)
