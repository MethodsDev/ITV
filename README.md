# Integrative Transcriptomics Viewer (ITV)

This Python package is based on a fork of the Genomeview package. It includes changes and helpers targeted at using ITV to view long reads sequencing data in the transcriptome context, and to allow more interactive comparisons.

See Tutorial/Vignette at [https://kinnex-documentation-external.readthedocs.io/en/latest/jupyter-notebooks/demo-ITV.html](https://kinnex-documentation-external.readthedocs.io/en/latest/jupyter-notebooks/demo-ITV.html)

## How to install

Requires conda

```
git clone https://github.com/MethodsDev/ITV.git && cd ITV
make install
```

You can then activate the conda environment with:
```
conda activate itv-1.1.2
```


## Running ITV in a Jupyter notebook

There are 2 options for running ITV in a Jupyter Notebook. You can either run Jupyter inside the ITV environment, or link the ITV environment to Jupyter installed in a different environment.

### Option 1: Run Jupyter server inside the ITV environment

```
conda activate itv-1.1.2 
jupyter lab
```

### Option 2: Run Jupyter server in a separate preexisting environment and link the ITV kernel to be usable

```
# link the ITV kernel
conda activate itv-1.1.2
ipython kernel install --user --name=itv-1.1.2 --display-name "Python ITV 1.1.2" # configure Jupyter to use Python kernel

# Then run jupyter from the system installation or a different conda environment:
conda activate my_other_env_with_jupyter          # this step can be omitted by using a different terminal window than before
jupyter lab
```

### Remote Jupyter session
If you are running Jupyter remotely (like on a VM), you can start it like this:
```
jupyter lab --no-browser --ip="0.0.0.0" --notebook-dir=/mnt/ --port=7777   # replace /mnt/ to base directory you want to use
```

You can then SSH to your VM and forward port 7777. Using a GCP VM, the command would look like this:
```
gcloud compute ssh MY_VM --project=PROJECT_NAME --zone=VM_ZONE --tunnel-through-iap -- -NL 7777:localhost:7777
```
