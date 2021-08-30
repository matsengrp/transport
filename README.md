# Comparing TCR repertoires via optimal transport

Branden J Olson, Stefan A Schattgen, Paul G Thomas, Philip Bradley, Frederick A Matsen IV

This repository contains the necessary dependencies and in-house code used to produce all analyses and figures as seen in the [Manuscript](https://github.com/matsengrp/transport-tex).
Below we describe how the installation of dependencies 
and provide instructions for one to reproduce all analyses. 

## Installation and Dependencies

1. First, clone the repository and it's sub modules

```bash
git clone --recurse-submodules https://github.com/matsengrp/transport.git
```

2. Next, you'll want to build the TCRDist executable

```bash
./build-pudtcrs.sh
```

3. Set up a proper python environment.

We recommend using [conda](https://docs.conda.io/en/latest/) for installation of python packages. You may use the spec-file.txt to create a mirror environment that currently works
```
conda create --name transport --file spec-file.txt
```

Other Requirements - running the code requires a few other dependencies called upon by our scripts

1. R >= 4.0.3
    additionally, you will need these third party packages for R:
    A. `install.package("estimatr")
    B. `install.package("RcmdrMisc")`

2. HMMER >= 3.2.1 #http://hmmer.org/
3. mafft v >= 7.453

On the Fred Hutch `rhino` servers, these modules are loaded like so:

```
module load R/4.0.3-foss-2020b
module load RepeatMasker/4.0.8-foss-2018b-Perl-5.28.0-HMMER
module load MAFFT/7.453-GCC-8.3.0-with-extensions
```


## Running manuscript analyses

The manuscript analyses are primarily contained within the [analyses directory](analyses/).
There, you will find a [README](analyses/README.md) that describes what each script does, as well as which script must be run first.
We provide an example of how one might run all analyses
with the order specified in [run\_all\_analyses.sh](run_all_analyses.sh).
The structure for which output is made
depends on [config.json](config.json).
Here, you may change where output get written to, as well as tweak common parameters
used for the manuscript. 

Once the analyses is run, you may explore the output directories to view various
raw data files and plots. 
For example, with the current configurations running
```
» bash run_all_analyses.sh
```
Taking roughly 30-40 minutes on a single core, this will produce a directory `output/` with the following structure:
```
output
├── cluster_iels
├── csv
├── dist_matrices
├── hmm
├── iel_clusters
├── json
├── mds
└── z_score
```
These output here is described in some detail along with the description of each script
in [analyses/README.\m\d](analyses/README.md).


### Code Description
The main code for the transport package lies in the `python` directory, which contains various modules.
The core of the package lies in the `TCRScorer` module, which takes two files `file_1` and `file_2` of TCR sequences as input, and computes the loneliness scores (called "enrichments" in this repo) of all TCRs in the repertoire corresponding to `file_2`.
This is then called by `TCRClusterer` to compute the loneliness cluster of TCRs, which in turn is called by `TCRMultiClusterer` to compute the top _k_ clusters for some specified _k_.
Thus, a generic analyses will suffice to call `TCRMultiClusterer` on two sequence files to obtain cluster inferences along with various files that include auxiliary information for each cluster.
There is also a `RandomizationTest` module which can be used to obtain significance estimates for the loneliness scores, although this has not yet been incorporated into `TCRScorer` or `TCRMultiClusterer`.

The `analyses` directory contains various scripts that were used to obtain the results in the transport manuscript.
These analyses can serve as illustrative examples of how the package can be used for new users.
For example, `analyses/combined_replicates.py` computes the loneliness scores used to obtain the OT-Tremont, OT-Revere, and OT-Ida clusters discussed in the manuscript.
These scripts usually write their results to a location within the root-level `output` directory.

The `R` directory contains scripts used to post-process the results from the `analyses` scripts, and generate plots.
These scripts also usually write their results to a location within `output`.
