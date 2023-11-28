# Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.1.0]

### Added

- support for input query formats: vcf
- --output-format option
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
