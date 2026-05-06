
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


## [1.4] - 2026-05-06

### Added
- Added support for 9 new species:
  - *Xenopus tropicalis* (UCB_Xtro_10.0)
  - *Gallus gallus* (bGalGal1.mat.broiler.GRCg7b)
  - *Sus scrofa* (Sscrofa11.1)
  - *Macaca mulatta* (Mmul_10)
  - *Caenorhabditis elegans* (WBcel235)
  - *Saccharomyces cerevisiae* (R64-1-1)
  - *Solanum lycopersicum* (SL3.0)
  - *Zea mays* (Zm-B73-REFERENCE-NAM-5.0)
  - *Glycine max* (Glycine_max_v2.1)
- Added `FIRST_CHR_CHECK` dict in `resources.py` mapping each species to its correct
  first chromosome identifier for existence checks (e.g. Roman numerals for yeast and
  *C. elegans*, `"2L"` for *Drosophila*).

### Changed
- Changed exon junction key separator from `"_"` to `"|"` in `ensembl.py`
  (`get_transcripts_dict`), `chooseTarget.py` (`format_junctions`, `get_junction_len`),
  and `construct_cdna.py`. This fixes junction key parsing for any species whose exon
  IDs contain underscores (e.g. yeast `YFL039C_mRNA-E1`, *C. elegans* `ZK617.1a.1.e1`).
- Updated chromosome sequence file reading in `construct_cdna.py` and `ensembl.py` to
  strip FASTA headers and join line-wrapped sequences into a single string before
  slicing. Ensembl FASTA files use 60-character line wrapping; passing wrapped sequences
  to Primer3 caused a `ValueError: Input line with no '='` crash.
- Replaced hardcoded `n = 1` / `n = "2L"` logic in `MASKED_SEQS()` with a
  `FIRST_CHR_CHECK` lookup, fixing spurious re-downloads on every run for species
  with non-numeric chromosome identifiers.

### Fixed
- Fixed `TypeError: cannot unpack non-iterable NoneType object` in `chooseTarget.py`
  for yeast and *C. elegans*, caused by `split("_")` shredding exon IDs that contain
  underscores (root cause of the junction key separator change above).
- Fixed `ValueError: Input line with no '='` from Primer3 for all non-vertebrate
  species whose chromosome `.txt` files are stored in FASTA format with line-wrapped
  sequences.
- Fixed `MASKED_SEQS()` triggering a re-download on every run for species using Roman
  numeral chromosome names (*S. cerevisiae*, *C. elegans*) because the existence check
  was always looking for `_1.txt` regardless of species.
