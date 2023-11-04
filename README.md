# interpolate-genetic-position

## Brief Summary

Use reference genetic map to interpolate genetic position for a query set of variants

## Overview

Certain statistical genetic applications require the analyst to annotate variants
by genetic position, typically in centimorgans or morgans. Unlike physical position,
which is canonically treated as a constant across genomes, genetic position may
potentially vary due to a number of more complicated factors. As such, it is typically
estimated at a nonuniform mesh of locations across the genome in a dedicated study,
and then those estimates are used to create estimates for individual downstream analyses.

In the case that the original estimation mesh overlaps with the experimental variant set,
annotation with genetic position is straightforward. However, often one must annotate
positions that were not present in the original mesh. To do so, canonically a parameter
in units (centimorgans / megabase) is reported across the mesh; this represents the rate
of change of genetic distance across physical distance. This is almost always a linear
approximation of a nonlinear function. To annotate new sites, the analyst must locate
the nearest mesh location with lower physical position than the experimental site,
and compute the genetic position for the site using the linear approximation


```
experimental.gpos = mesh.gpos + (experimental.ppos - mesh.ppos) / 1000000 * mesh.rate
```

where `gpos` is genetic position in centimorgans, `ppos` is physical position in the
reference genome, and `rate` is the (centimorgans / megabase) parameter mentioned above.

This calculation is fairly straightforward, though minor complications arise in corner cases:

* an experimental variant physically located *before* the first mesh position on a chromosome
  is assigned the genetic position 0; this is true of all such variants meeting this condition
* similarly, an experimental variant physically located *after* the last mesh position on
  a chromosome is assigned the genetic distance of that last mesh position
  * since the first mesh location is arbitrarily assigned genetic position 0, these are actually
    the same condition, but they tend to look different at first glance
* experimental genetic position is typically only computed for autosomes and the X chromosome
* experimental genetic position is genome build-specific; there are conceptually different
  estimates depending on which reference genome is used in the estimate
  * note, however, that many versions of the human genetic recombination map were
    simply translated from earlier genome reference builds with [liftOver](https://genome.ucsc.edu/cgi-bin/hgLiftOver),
    which has substantial flaws

All of the above is manageable, but fiddly and error prone, and is a repeated task that
people implement individually each time they do it.

This repo contains a standalone utility that performs this annotation and interpolation
on a variety of different input data formats and using an assortment of possible genetic maps.
It is (will be) tested such that its output is reliable, and when approximations are required,
those will be described and exposed to user configuration as much as possible.



## Installation

The recommended installation method for this package _will be_ via my conda channel; but while
this is in pre-release development, interested users will have to build from source.

### Requirements

See the provided [environment.yaml](environment.yaml) file for conda-formatted dependencies for
building this package.

### Build

By default, a build process involving a [conda/mamba](https://mamba.readthedocs.io/en/latest/installation.html) environment is supported.

  - if you wish to use `conda` and it's not currently available, you can install it with the instructions [here](https://mamba.readthedocs.io/en/latest/mamba-installation.html#mamba-install)
  - navigate into your project directory (interpolate-genetic-position)
  - create the `conda` environment for installation as follows:
  
     `mamba env create -f environment.yaml`
  - activate the conda environment:
  
     `mamba activate interpolate-genetic-position-env`
  - update (create) the necessary `configure` scripts with `autoreconf`:
  
     `autoreconf --force --install`
	 
     - note that this can also be run with `./generate.bash` inside the repo
  - run `configure`:
  
	 `./configure --with-boost=${CONDA_PREFIX} --with-boost-libdir=${CONDA_PREFIX}/lib`

	 - if you are planning on installing software to a local directory, run instead `./configure --prefix=/install/dir [...]`
  - run `make -j{ncores}`

  - if desired, run `make install`. if permissions issues are reported, see above for reconfiguring with `./configure --prefix`.
  
## Usage

By default, the final compiled program can be run with

`./interpolate-genetic-position.out`


### Command Line Parameters

|Parameter|Description|
|---|---|
|`--input`<br>`-i`|Input file of variants or regions to annotate. Needs to be sorted, chromosome and position. Can be gzipped. If not specified, will be read from stdin.|
|`--preset`<br>`-p`|Format of input variant file. Accepted formats: `bim`, `map`, `bed`.|
|`--genetic-map`<br>`-g`|Input recombination map. Needs to be sorted, chromosome and position. Can be gzipped.|
|`--map-format`<br>`-m`|Format of recombination map. Accepted formats: `bolt`, `bedgraph` (see below for further discussion).|
|`--output`<br>`-o`|Output file. Will match format of input. Cannot currently be gzipped. If not specified, will be written to stdout.|
|`--verbose`<br>`-v`|Whether to print extremely verbose debug logs. You probably don't want this.|


## How to Choose a Recombination Rate File

There are an assortment of different recombination rate/genetic map files available publicly online.
The choice of which map to use has a variety of constraints:

* Consider your input experimental data. Depending on the genome build used, you may be constrained
  in which map you choose. The oldest of these maps dates back to HapMap2/hg18. A large number of maps
  in later builds are derived directly from this map with liftOver, and are arguably of low quality
* Certain experimental applications might benefit from the use of a recombination map that was estimated
  from a population that matches the genetic ancestry of your dataset. In general, there haven't been
  many great findings derived from using population-specific maps, and most people tend to rely on
  either cosmopolitan or single ancestry (read: European) maps. Most tools don't really seem to care too much

The following primary maps are supported by this tool:

* Genetic maps released with BOLT-LMM (`--map-format bolt`)
  * Available in the tarball [in this directory](https://alkesgroup.broadinstitute.org/BOLT-LMM/downloads/). Download
    and extract the latest version (2.4.1 as of 2023); the map you want is in the `tables` subdirectory.
    * Make sure you choose the one matching your dataset's genome build
  * These maps cover hg17, hg18, hg19 (GRCh37), and hg38 (GRCh38).
* Recombination rate tracks from UCSC (`--map-format bedgraph`)
  * Available as bigwigs [in this directory](https://hgdownload.soe.ucsc.edu/gbdb/hg38/recombRate/). Download
    the version you prefer (I have used [this one](https://hgdownload.soe.ucsc.edu/gbdb/hg38/recombRate/recombAvg.bw)).
  * The link in this readme is just for hg38/GRCh38. You can go digging around for others if you like; they all
    have benefits and drawbacks. Hopefully the track format is consistent enough that it'll be compatible.
  * This program does not currently have the ability to directly read bigwigs, and furthermore needs the chromosomes
    in the map to be sorted. To use the UCSC track with this program, do some version of the following:
    * `mamba install -c bioconda -c conda-forge ucsc-bigwigtobedgraph`
    * `bigWigToBedGraph recombAvg.bw recombAvg.bedgraph`
    * `sed 's/^chrX/23/ ; s/^chr//' recombAvg.bedgraph | sort -k 1,1g -k 2,2g -k3,3g | sed -r 's/^23/chrX/ ; s/^([^c])/chr\1/' > recombAvg.sorted.bedgraph`



## Version History

28 10 2023: project generated from cookiecutter template
