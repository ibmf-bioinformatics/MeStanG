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
You can install these packages using a package manager like micromamba as follows:

```bash=1
micromamba create -n metasim python=3.8 biopython=1.83 numpy=1.26.4 -c conda-forge
```

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

For detailed instructions on usage see the [Usage Manual](manuals/Manual.md).

## Credits

MeStanG is developed in Andres S. Espindola lab at the Institute for Biosecurity and Microbial Forensics - Oklahoma State University

Main code contributors:

* MeStanG modules and Manual: Daniel Ramos Lopez

## Cite

Use the following doi to cite the usage of this tool: https://doi.org/10.3390/biology14010069

## Issues and help

Users can send any questions about MeStanG usage on the issues tab, before submitting an issue please consider looking through the manual and browsing existing issues.

## Acknowledgements

PhD. Andres S. Espindola - Oklahoma State University

PhD. Francisco Flores Flor - Universidad de las Fuerzas Armadas - ESPE

Contributors to the development and users of this tool
