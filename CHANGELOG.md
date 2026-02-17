
# Change Log
All notable changes to this project will be documented in this file.
 
The format is based on [Keep a Changelog](http://keepachangelog.com/)
and this project adheres to [Semantic Versioning](http://semver.org/).
 
## [Unreleased] - 2022-12-16
 
Here we write upgrading notes for brands. It's a team effort to make them as
straightforward as possible.
 
### Added
- Added the setup.py file
- Created new function to download the files requeriments
 
### Changed
 - Moved the pipe to exon_surfer.py
 
### Fixed
 - Fixed the package structure

 ## [0.1.0] - 2021-12-21

### Added
- Added the "create_index_table" function, to match the transcript ID with the
  gene symbol

### Changed

### Fixed
- Fixed the "make_blast_db" function, to create the blast database with the
  correct cDNA file

## [1.0] - 2023-12-18

### Fixed
- Fixed the if else to obtain the junctions in CreatePrimers


## [1.3] - 2026-02-17

### Added
- Added support for *Litomosoides sigmodontis* (WormBase ParaSite) as a custom genome.
- Added `get_genome_data` function in `resources.py` to centralize genome loading strategies.
- Added Zenodo resource links for *L. sigmodontis* genomic and blast databases.

### Changed
- Refactored `CreatePrimers` to dynamically load genome data via `resources.py` instead of hardcoded Ensembl calls.
- Updated dimer filtering logic: the filter is now skipped if it eliminates all primer candidates, preventing empty result crashes.

### Fixed
- Fixed "unhashable type: list" error in `chooseTarget.py` for species without defined canonical transcripts (e.g., *Arabidopsis thaliana*).
- Fixed "Expected a 1D array" error in specific genes where the dimer filter resulted in an empty dataframe.
