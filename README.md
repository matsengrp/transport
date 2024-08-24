# transport

This repository accompanies the manuscript 
[_Comparing TCR repertoires via optimal transport_](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1010681), 
by Olson, Schattgen, Thomas, Bradley, and Matsen.
It provides instructions to reproduce the analyses.

Note that this repository contains data generously shared by our collaborators
and described in the manuscript _Intestinal Intraepithelial Lymphocyte
Repertoires are Imprinted Clonal Structures Selected for MHC Reactivity_ by
Schattgen, Crawford, Van de Velde, Chu, Mazmanian, Bradley, and Thomas.
Please cite that paper if you use these data.

**This package is no longer under active development.**
The [IRTransport](https://github.com/zacmon/ir_transport) package provides an 
updated, convenient, more efficient, and flexible recoding of the `transport` 
package.

## Installation and Dependencies

1. First, clone the repository and its sub modules

```bash
git clone --recurse-submodules https://github.com/matsengrp/transport.git
```

2. Next, you'll want to build the TCRDist executable

```bash
cd transport
./build-pubtcrs.sh
```

3. Set up a proper Python environment

We recommend using [conda](https://docs.conda.io/en/latest/) for installation of python packages:
```bash
# will name the conda env 'transport'
conda env create --file environment.yml
```

4. Install other requirements

    * R >= 4.0.3 with additional packages:
        `install.packages(c("estimatr", "RcmdrMisc", "rjson", "segmented"))`
    * HMMER >= 3.2.1
    * mafft >= 7.453

On the Fred Hutch `rhino` servers, these modules are loaded like so:

```bash
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
bash run_all_analyses.sh
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
in [analyses/README.md](analyses/README.md).


### Code Description
The main code for the transport package lies in the `python` directory, which contains various modules.
The core of the package lies in the `TCRScorer` module, which takes two files `file_1` and `file_2` of TCR sequences as input, and computes the loneliness scores (called "enrichments" in this repo) of all TCRs in the repertoire corresponding to `file_2`.
There is also a `RandomizationTest` module which can be used to obtain significance estimates for the loneliness scores, although this has not yet been incorporated into `TCRScorer` or `TCRMultiClusterer`.

The `analyses` directory contains various scripts that were used to obtain the results in the transport manuscript.
These analyses can serve as illustrative examples of how the package can be used for new users.
For example, `analyses/combined_replicates.py` computes the loneliness scores used to obtain the OT-Tremont, OT-Revere, and OT-Ida clusters discussed in the manuscript.
These scripts usually write their results to a location within the root-level `output` directory.

The `R` directory contains scripts used to post-process the results from the `analyses` scripts, and generate plots.
These scripts also usually write their results to a location within `output`.
