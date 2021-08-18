# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [2.0.0] - 2021-08-18
### Added
  - Method for scoring modules using a DGEList object from edgeR

### Changed
  - Simplified methods
  - Added package namespace declarations to all functions
  - Replace {magrittr} pipe with native R pipe where possible
  - Remove all instances of `%<>%`
  - Remove all instances of `return()`
  - `prepGMT` now looks for the `term` column instead of `ont`

[2.0.0]: https://github.com/milescsmith/moduleScoreR/releases/tag/2.0.0
