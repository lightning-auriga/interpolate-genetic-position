# Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.2.0]

### Changed

- in experimental mode, when genetic position is incremented between bed query regions
  by a fixed genetic distance, the program now refrains from incrementing between
  consecutive query regions with the same bedfile column 4. conceptually, this allows
  multiple distinct regions to be assigned to the same conceptual unit for the purposes
  of this experimental mode. the applications of this functionality are very niche. however,
  when the increment is set to 0, as should be the case for most users, this change
  will have no impact whatsoever (see #16)

### Fixed

- bedfile queries with non-contiguous end/start positions on the same chromosome
  cause an additional placeholder entry to be injected into the output corresponding
  to changes in the input genetic map's rate. this prevents incorrect interpolation
  of genetic position from the output genetic map when the uncovered region is
  implicitly assumed to have the previous map region's rate (see #15)

## [1.1.0]

### Added

- support for input query formats: vcf
- --output-format option
- output "bed" format officially changed to "bolt"
- generic (as opposed to C++ specific) linters
- htslib build dependency
- improved test coverage

### Changed

- input and output format are no longer required to be the same
  - there are still restrictions on sane conversions

### Fixed

- verbose logging no longer interferes with output streaming


## [1.0.0]

### Added

- support for input query formats: bed, bim, map, snp
- support for genetic maps: bolt, bedgraph, bigwig
- support for output in morgans
- experimental support for genetic position increments at boundaries of bed query regions

[//]: # [Unreleased]

[//]: # (- Added)
[//]: # (- Changed)
[//]: # (- Deprecated)
[//]: # (- Removed)
[//]: # (- Fixed)
[//]: # (- Security)
