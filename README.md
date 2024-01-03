# interpolate-genetic-position


[![CircleCI](https://dl.circleci.com/status-badge/img/gh/lightning-auriga/interpolate-genetic-position/tree/default.svg?style=svg)](https://dl.circleci.com/status-badge/redirect/gh/lightning-auriga/interpolate-genetic-position/tree/default)

[![codecov](https://codecov.io/gh/lightning-auriga/interpolate-genetic-position/graph/badge.svg?token=l5pw6XfJ7l)](https://codecov.io/gh/lightning-auriga/interpolate-genetic-position)

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

There are two primary installation methods, depending on your needs.

### With conda (recommended)

  - If needed, [install mamba](https://mamba.readthedocs.io/en/latest/mamba-installation.html#mamba-install)
  - Install the package with mamba:

    `mamba create -n igp -c https://raw.githubusercontent.com/lightning-auriga/conda-builds/default/conda-builds -c bioconda -c conda-forge interpolate-genetic-position`
  - Activate the resulting environment:

    `mamba activate igp`
  - The tool should now be available as `interpolate-genetic-position.out`

### Manual Build

#### Requirements

See the provided [environment.yaml](environment.yaml) file for conda-formatted dependencies for
building this package.

#### Build

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

`interpolate-genetic-position.out`


### Command Line Parameters

|Parameter|Description|
|---|---|
|`--input`<br>`-i`|Input file of variants or regions to annotate. Needs to be sorted, chromosome and position. Can be gzipped. If not specified, will be read as plaintext from stdin.|
|`--preset`<br>`-p`|Format of input variant file. Accepted formats: `bim`, `map`, `snp`, `vcf`, `bed`.|
|`--genetic-map`<br>`-g`|Input recombination map. Needs to be sorted, chromosome and position. Can be gzipped (except bigwigs). If not specified, will be read as plaintext from stdin.|
|`--map-format`<br>`-m`|Format of recombination map. Accepted formats: `bolt`, `bedgraph`, `bigwig` (see below for further discussion).|
|`--output`<br>`-o`|Output file. Will match format of input. Cannot currently be gzipped. If not specified, will be written to stdout.|
|`--output-format`<br>`-f`|Format of output file. Accepted formats: `bolt`, `bim`, `map`, `snp (see below for further discussion).|
|`--verbose`<br>`-v`|Whether to print extremely verbose debug logs. You probably don't want this.|
|`--output-morgans`|Report output genetic position in morgans, instead of the default centimorgans.|
|`--region-step-interval`|Add a fixed genetic distance at the boundaries of end positions of bedfile region queries, such that the output data have a step-like structure. This functionality is included for experimental purposes, and in most applications this setting should be kept at its default of 0.|
|`--help`<br>`-h`|Print brief help message and exit.|
|`--version`|Print version string for current build.|

Note that, of the above, either `-i` or `-g` can be read from stdin, but not both.

## Valid Combinations of Input and Output Formats

This program can attempt to automatically reformat input files into different format output
files, but the effectiveness of such a conversion with the available information varies.

|Input Format|Valid Output Formats|Notes|
|---|---|---|
|bim|bim, map, snp||
|map|map|Map files lack allele information, and so allele-containing formats are not possible.|
|snp|bim, map, snp||
|vcf|bim, map, snp|Vcf is not itself a supported output format. Note that for markers with multiple alternate alleles, only the first will be reported.|
|bed|bolt|Input bed regions are converted into bolt-format genetic maps.|



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
* Recombination rate tracks from UCSC (`--map-format bigwig` or `--map-format bedgraph`)
  * Available as bigwigs [in this directory](https://hgdownload.soe.ucsc.edu/gbdb/hg38/recombRate/). Download
    the version you prefer (I have used [this one](https://hgdownload.soe.ucsc.edu/gbdb/hg38/recombRate/recombAvg.bw)).
  * The link in this readme is just for hg38/GRCh38. You can go digging around for others if you like; they all
    have benefits and drawbacks. Hopefully the track format is consistent enough that it'll be compatible.
  * In theory, this should work on any of the aforementioned genome builds; in practice, I've only tested it
    for hg38 so far.
  * If you have mismatched contig issues with a bigwig, you can convert a bigwig to bedgraph and have it in any chromosome annotation.
    However, the bedgraph needs to be sorted by karyotyping order, not the lexicographical sort the UCSC tracks use.
    To convert a UCSC track to a bedgraph for use with this program, do some version of the following:
    * `mamba install -c bioconda -c conda-forge ucsc-bigwigtobedgraph`
    * `bigWigToBedGraph recombAvg.bw recombAvg.bedgraph`
    * `sed 's/^chrX/23/ ; s/^chr//' recombAvg.bedgraph | sort -k 1,1g -k 2,2g -k3,3g | sed -r 's/^23/chrX/ ; s/^([^c])/chr\1/' > recombAvg.sorted.bedgraph`
  * Note that there's one potentially significant advantage to using bedgraphs: as currently implemented,
    the program will load all recombination data for a chromosome at once (clearing memory between chromosomes).
    For most recombination maps and computers, this RAM usage is pretty trivial. However, if you have a particularly
    fine-grained recombination map, the intervals themselves could be a memory block. If you convert such a map to
    bedgraph (and sort it) in advance, the program will stream linewise and have fixed, constant, low memory use.

## Notes Specific to Interval (bedfile) Input

For idiosyncratic reasons, input queries in bed format are emitted as bolt-format genetic maps. This is due to the need
to specify both rate and position in the output blocks. Note that query blocks are fragmented based on whatever
mosaic pattern in which they happen to intersect with the initial genetic map. Due to the lack of end positions in
bolt format maps, a final bolt entry per-chromosome is injected into the output with fixed rate 0, to prevent implied
interpolation with the rate of the last estimated block.

The flag `--region-step-interval [double-precision value]` causes the bolt-format genetic map emitted for
bedfile queries to contain stepwise increments in genetic position at the end position of each query interval.
This is for experimental purposes; for most practical uses of this tool, this flag should be left at its default of 0.

When `--region-step-interval` is non-zero, bedfile column 4 is considered a group label, such that consecutive queries
with the same column 4 entry will not have their genetic position incremented by the fixed specified interval.
This behavior allows distinct query regions to be assigned to the same conceptual unit. This behavior is obviously very niche,
and since `--region-step-interval` should almost always be 0, this behavior will have no practical impact in most circumstances.

## I/O Streams

This program can accept input files as streams and can emit output to stream. To have a file be read from
or written to stream, simply omit its relevant flag (`-i`, `-g`, `-o`). For input streams, the corresponding
format flag should still be set. Only one input stream (either `-i` or `-g`) can be read from an input stream
per run.

## Example Use Cases

### Take a plink-format genotype dataset, annotate with interpolated genetic positions

This is the most basic expected use case.

```bash
interpolate-genetic-position.out -i infile.bim -p bim -g genetic_map.tsv -m bolt -o infile_with_gpos.bim -f bim
```

### Take a vcf, annotate markers with interpolated genetic positions

Note that vcf **output** format is not supported, so you'll have to do some further processing if you want to add
the interpolated values back into e.g. the INFO field somewhere.

```bash
interpolate-genetic-position.out -i infile.vcf.gz -p vcf -g genetic_map.tsv -m bolt -o interpolated_values.map -f map
```

### Take a bed file of regions, turn it into its own recombination map

Regions don't have a single genetic distance associated with them, and so for completeness, when annotating
a bed region file, the output becomes a genetic map with associated genetic distance and rate (cM/Mb) data.

```bash
interpolate-genetic-position.out -i infile.bed -p bed -g genetic_map.tsv -m bolt -o new_recombination_map.tsv -f bolt
```

### Take a bed file of regions, turn it into its own recombination map, and use that map to annotate a second file

Either the input query file (`-i`), or the input recombination map (`-g`), but not both at once, can be read
from stream by leaving the corresponding argument unspecified. Note that the corresponding format flag
is still required. The output can also be streamed in the same manner with the same restriction.

```bash
interpolate-genetic-position.out -i infile.bed -p bed -g genetic_map.tsv -m bolt -f bolt |
interpolate-genetic-position.out -i infile.vcf -p vcf -m bolt -o annotated_variants.map -f map
```

## Version History

See ChangeLog.md.
