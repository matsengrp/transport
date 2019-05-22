# Comparing TCR repertoires via optimal transport

## Installation
We recommend using [conda](https://docs.conda.io/en/latest/) for installation.
It should be possible to get things running using other package managers like pip via analogous commands minus the conda-specific environment setup, but we haven't tested this.

Running the following commands will set up the conda environment with the necessary packages:
```
conda update -y conda
conda create -n "transport" python=3
activate transport
conda install -y matplotlib numpy pandas seaborn
conda install -y -c conda-forge pot
```

Then, running `jupyter notebook` should start a tab in your browser.

For stoat, running
```
jupyter notebook --no-browser --port=8889 --ip=0.0.0.0
```
and copying the URL that is output into your browser (e.g. `http://stoat.8891/?token=LONG_STRING_OF_NUMBERS`) should do the trick if you're conncted to Marconi.
(I haven't figured out how to interact with a jupyter notebook remotely through the VPN.)
