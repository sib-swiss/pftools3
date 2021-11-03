[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/pftools/README.html)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/pftools/badges/license.svg)](https://anaconda.org/bioconda/pftools)
[![Anaconda-Server Badge](https://img.shields.io/conda/dn/bioconda/pftools.svg?style=flat)](https://anaconda.org/bioconda/pftools)

PfTools
=========================================

## Table of Contents

   * [Foreword](#foreword)
   * [Installation](#installation)
     * [Using Docker](#using-docker)
     * [Using Singularity](#using-singularity)
     * [Bioconda](#bioconda)
     * [Manually](#manually) 
   * [Generalized profile syntax](#generalized-profile-syntax)
   * [Algorithms description](#algorithms-description)
   * [Applications of the Pftools](#applications-of-the-pftools)
   * [Authors](#authors)

# Foreword

(C) Copyright SIB Swiss Institute of Bioinformatics
available from  https://github.com/sib-swiss/pftools3 under GPL v2. See LICENSE.


Version 3 contains the original FORTRAN 77 pftools (release 2.3)
and the new pftoolsV3 programs.

# Installation

### Using Docker

First you must have [Docker](https://docs.docker.com/get-docker/) installed and running.  
Secondly have a look at the availabe pftools biocontainers at [quay.io](https://quay.io/repository/biocontainers/pftools?tab=tags) or at [Docker Hub](https://hub.docker.com/r/sibswiss/pftools).  
Then:
  ```
# get the chosen pftools container version
docker pull quay.io/biocontainers/pftools:2.3.5--h4333106_0
   or
docker pull sibswiss/pftools:3.2.8
# use an pftools's tool e.g. pfscan 
docker run quay.io/biocontainers/pftools:2.3.5--h4333106_0 pfscan -h
   or
docker run sibswiss/pftools:3.2.8 pfscan -h
  ```

### Using Singularity

First you must have [Singularity](https://sylabs.io/guides/master/user-guide/quick_start.html) installed and running.
Secondly have a look at the availabe pftools biocontainers at [quay.io](https://quay.io/repository/biocontainers/pftools?tab=tags) or at [Docker Hub](https://hub.docker.com/r/sibswiss/pftools).  
Then:
```
# get the chosen pftools container version
singularity pull docker://quay.io/biocontainers/quay.io/biocontainers/pftools:2.3.5--h4333106_0
   or
singularity pull docker://sibswiss/pftools:3.2.8
# run the container
singularity run pftools_2.3.5--h4333106_0.sif
```

You are now in the container. You can use an pftools's tool e.g. pfscan doing 
```
pfscan -h
```

## Bioconda

```
conda install -c bioconda pftools
```

## Manually

See [here](./INSTALL) for more information

**Prerequisite**  

  * cmake >= 3.7
  * gcc   >= 4.6
  * gfortran (or g77 or f77)  for release 2.3 source code
  * perl  >= 5.5.3   
    * File::Slurp

**BUILD**  
```
git clone https://github.com/sib-swiss/pftools3.git
cd pftools3
mkdir build
cd build/
cmake ..
```

**COMPILE**  
```
make
```

**INSTALL**  
```
make install
```

**CHECK/TEST**  
```
make test
```

After installation, in the share/examples/ subdirectory, the *test_V3.sh* shell script is a good starting point for using pfsearchV3/pfscanV3.

# Generalized profile syntax

A description of the generalized profile syntax is given in file:

- [doc/profile.txt](https://raw.githubusercontent.com/sib-swiss/pftools3/master/doc/profile.txt)  (original document)
- [doc/profile.pdf](https://raw.githubusercontent.com/sib-swiss/pftools3/master/doc/profile.pdf)  (revised and completed version)

it was originally published in

* Bucher P, Bairoch A.
  A generalized profile syntax for biomolecular sequence motifs
  and its function in automatic sequence interpretation.
  Proc Int Conf Intell Syst Mol Biol. 1994;2:53-61.
  PubMed PMID: [7584418](https://www.ncbi.nlm.nih.gov/pubmed/7584418).


# Algorithms description

Technical details about how profiles can be constructed and parametrized
are summarized in file:

- [doc/profile.pdf](https://raw.githubusercontent.com/sib-swiss/pftools3/master/doc/profile.pdf)

The very first paper describing the PFTOOLS algorithms is

* Lüthy R, Xenarios I, Bucher P.
  Improving the sensitivity of the sequence profile method.
  Protein Sci. 1994 Jan;3(1):139-46.
  PubMed PMID: [7511453](https://www.ncbi.nlm.nih.gov/pubmed/7511453); PubMed Central PMCID: PMC2142471.

The generalized profile alignment method is closely related to other "classical"
algorithm for aligning sequences. For example, it encompasses the Smith-Waterman
algorithm and the Viterbi decoding of profile-HMM (as implemented in HMMER2 for
example). Relationships between these algorithm were investigated in

* Bucher P, Hofmann K.
  A sequence similarity search algorithm based on a probabilistic interpretation of an alignment scoring system.
  Proc Int Conf Intell Syst Mol Biol. 1996;4:44-51. Review.
  PubMed PMID: [8877503](https://www.ncbi.nlm.nih.gov/pubmed/8877503).

* Bucher P, Karplus K, Moeri N, Hofmann K.
  A flexible motif search technique based on generalized profiles.
  Comput Chem. 1996 Mar;20(1):3-23.
  PubMed PMID: [8867839](https://www.ncbi.nlm.nih.gov/pubmed/8867839).

Relatively detailed explanations about the profile normalized scores, as well as its
comparisons with other popular statistics for sequence alignments can be found in

* Pagni M, Jongeneel CV.
  Making sense of score statistics for sequence alignments.
  Brief Bioinform. 2001 Mar;2(1):51-67.
  PubMed PMID: [11465063](https://www.ncbi.nlm.nih.gov/pubmed/11465053).

The heuristic score is succinctly described in

* Schuepbach T, Pagni M, Bridge A, Bougueleret L, Xenarios I, Cerutti L.
  pfsearchV3: a code acceleration and heuristic to search PROSITE profiles.
  Bioinformatics. 2013 May 1;29(9):1215-7. doi: 10.1093/bioinformatics/btt129.
  PubMed PMID: [23505298](https://www.ncbi.nlm.nih.gov/pubmed/23505298); PubMed Central PMCID: PMC3634184.

# Applications of the Pftools

Two databases were created based on the PFTOOLS technology: PROSITE and HAMAP
and they are still actively maintained

1. https://prosite.expasy.org/
1. https://hamap.expasy.org/

The PFTOOLS were initially designed with handling capabilities of DNA sequences.
The latest released pfsearchV3 feature support for FASTQ and SAM formats. DNA
applications are for example given in

* Pagni M, Niculita-Hirzel H, Pellissier L, Dubuis A, Xenarios I, Guisan A, Sanders IR, Goudet J, Guex N.
  Density-based hierarchical clustering of pyro-sequences on a large scale - the case of fungal ITS1.
  Bioinformatics. 2013 May 15;29(10):1268-74. doi: 10.1093/bioinformatics/btt149.
  PubMed PMID: [23539304](https://www.ncbi.nlm.nih.gov/pubmed/23539304) ; PubMed Central PMCID: PMC3654712.

* Schmid-Siegert E, Richard S, Luraschi A, Mühlethaler K, Pagni M, Hauser PM.
  Mechanisms of Surface Antigenic Variation in the Human Pathogenic Fungus Pneumocystis jirovecii.
  MBio. 2017 Nov 7;8(6). pii: e01470-17. doi: 10.1128/mBio.01470-17.
  PubMed PMID: [29114024](https://www.ncbi.nlm.nih.gov/pubmed/29114024); PubMed Central PMCID: PMC5676039.

# Authors

Mas:
- Philipp Bucher developped the Fortran code
- Thierry Schuepbach developped the C code

Other contributors:
- Kay Hofmann
- Volker Flegel
- Edouard de Castro
- Lorenzo Cerruti
- Marco Pagni
- Sébastien Moretti
- Jerven Tjalling Bolleman

[SIB Swiss Institute of Bioinformatics](https://www.sib.swiss/)
[Vital-IT Group](https://www.vital-it.ch/)
Quartier Sorge - Batiment Amphipole
1015 Lausanne
Switzerland

