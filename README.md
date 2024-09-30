# MeStanG

**Me**tagenomic **Stan**dards **G**enerator (MeStanG) for HTS Nanopore datasets

Resource for simulating *de novo* nanopore datasets resembling systematic sampling sequencing data. 

The pipeline has been tested with Python 3.12.3, using the modules `biopython = 1.83` and `numpy = 1.26.4`. Tested on CentOS Linux 7.9.2009, though it is likely to run on any platform compatible with Python 3. 

## Installation

Download the most recent source code, uncompress it, and add all `.py` files to `PATH` and/or `PYTHONPATH`. You can also place the `.py` files on your `bin` folder.

## Requirements

Python 3.12.3 and the following modules:

```
biopython==1.83
numpy==1.26.4
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

## Issues and help

Users can send any questions about MeStanG usage on the issues tab, before submitting an issue please consider looking through the manual and browsing existing issues.

## Acknowledgements

PhD. Andres S. Espindola - Oklahoma State University

PhD. Francisco Flores Flor - Universidad de las Fuerzas Armadas - ESPE

Contributors to the development and users of this tool
