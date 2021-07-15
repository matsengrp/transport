# Comparing TCR repertoires via optimal transport

## Environment and Dependencies

Clone the repository and it's submodules

using ssh keys
```bash
git clone --recurse-submodules git@github.com:matsengrp/transport.git
```
using https
```bash
git clone --recurse-submodules https://github.com/matsengrp/transport.git
```

First, you'll want to build the TCRDist executable
```bash
./build-pudtcrs.sh
```

Next, set up a proper python environment.
We recommend using [conda](https://docs.conda.io/en/latest/) for installation of python packages. You may use the spec-file.txt to create a mirror environment that currently works

```
conda create --name transport --file spec-file.txt
```

Other Requirements - running the code requires a few other dependencies called upon by our scripts

1. R >= 4.0.3
2. HMMER >= 3.2.1 #http://hmmer.org/
3. mafft v >= 7.453

on the Fred Hutch `rhino` servers, I load these modules like so:
```
module load R/4.0.3-foss-2020b
module load RepeatMasker/4.0.8-foss-2018b-Perl-5.28.0-HMMER
module load MAFFT/7.453-GCC-8.3.0-with-extensions
```


## Running manuscript analysis

Next, we'll walk you through running the three primary manuscript analysis and their respective input/outputs 

### IEL's data

.... TODO

### Code Description
The main code for the transport package lies in the `python` directory, which contains various modules.
The core of the package lies in the `TCRScorer` module, which takes two files `file_1` and `file_2` of TCR sequences as input, and computes the loneliness scores (called "enrichments" in this repo) of all TCRs in the repertoire corresponding to `file_2`.
This is then called by `TCRClusterer` to compute the lonelinest cluster of TCRs, which in turn is called by `TCRMultiClusterer` to compute the top _k_ clusters for some specified _k_.
Thus, a generic analysis will suffice to call `TCRMultiClusterer` on two sequence files to obtain cluster inferences along with various files that include auxiliary information for each cluster.
There is also a `RandomizationTest` module which can be used to obtain significance estimates for the loneliness scores, although this has not yet been incorporated into `TCRScorer` or `TCRMultiClusterer`.

The `analyses` directory contains various scripts that were used to obtain the results in the transport manuscript.
These analyses can serve as illustrative examples of how the package can be used for new users.
For example, `analyses/combined_replicates.py` computes the loneliness scores used to obtain the OT-Tremont, OT-Revere, and OT-Ida clusters discussed in the manuscript.
These scripts usually write their results to a location within the root-level `output` directory.

The `R` directory contains scripts used to post-process the results from the `analyses` scripts, and generate plots.
These scripts also usually write their results to a location within `output`.
