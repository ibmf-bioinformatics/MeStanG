# MeStanG

**Me**tagenomic **Stan**dards **G**enerator (MeStanG) for HTS Nanopore datasets

Resource for simulating *de novo* nanopore datasets resembling systematic sampling sequencing data. 

The pipeline has been tested with Python 3.12.3, using the modules `biopython = 1.83` and `numpy = 1.26.4`. It was tested on CentOS Linux 7.9.2009, though it will likely run on any platform compatible with Python 3. Multithreading managed by multiprocessing on Python is compatible with Linux but it might not work properly on Windows and MacOS, single-thread runs are compatible with all OS.

## Installation

Add the directory containing `MeStanG.py` to `PATH` and/or `PYTHONPATH` or place `MeStanG.py` in your `bin` folder. Make sure the `.py` is executable.

## Requirements

Python 3.12.3 and the following modules:

```
biopython==1.83
numpy==1.26.4
```
You can install these packages using a package manager like [micromamba](https://github.com/mamba-org/micromamba-releases) as follows:

```bash=1
micromamba create -n mestang python=3.12.3 biopython=1.83 numpy=1.26.4 -c conda-forge
```

Now, to activate the environment:

```bash=1
micromamba activate mestang
```

Now you can use MeStanG by either adding it to your PATH, or a workaround is to copy it to your current environment `bin` folder so every time you load the environment, you can use it. First, you need to know where your environment is stored, so you can run

```bash=1
micromamba env list
ls <environmentlocation>/bin
cp MeStanG.py <environmentlocation>/bin
chmod +x <environmentlocation>/bin/MeStanG.py
```
Now you can use MeStanG.py

## Usage

```Python=22
usage: MetaStaG.py [-h] [-v] {env,host} ...

        MeStanG
        -----------------------------------------------------------
        Given a set of genomes, assemblies, or contigs
        generate standard datasets of raw ONT reads

options:
  -h, --help     show this help message and exit
  -v, --version  show program's version number and exit

subcommands:

  For detailed usage of each sample source:
      MeStanG.py sample -h
  -------------------------------------------------------

  {env,host}     You may run on env or host sample source
    env          Environmental samples
    host         Host/pathogen samples
```

For detailed instructions on usage, see the [Usage Manual](manuals/Manual.md).

## Credits

MeStanG is developed in the Andres S. Espindola lab at the Institute for Biosecurity and Microbial Forensics - Oklahoma State University

Main code contributors:

* MeStanG modules and Manual: Daniel Ramos Lopez

## Cite

Use the following doi to cite the usage of this tool: https://doi.org/10.3390/biology14010069

## Issues and help

Users can send any questions about MeStanG usage to the issues tab. Before submitting an issue, please consider reviewing the manual and browsing existing issues.

## Acknowledgements

PhD. Andres S. Espindola - Oklahoma State University

PhD. Francisco Flores Flor - Universidad de las Fuerzas Armadas - ESPE

Contributors to the development and users of this tool
