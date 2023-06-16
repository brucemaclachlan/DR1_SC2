# Structural definition of HLA class-II presented SARS-CoV-2 epitopes reveals a mechanism to escape pre-existing CD4+ T cell immunity 
This is a repository that relates to the paper Chen et al 2023 in Cell Reports

## Description

Using this repository, all structural figures from the above manuscript can be recreated from coordinate (.pdb) and structure factor (.mtz) files available via the PDB.

The analysis consists of two scripts
* pMHC_analyse_v2.py

This creates the general overviews of each structure as well as performs omit map analysis
* pocket_analysis_v2.py

This creates the views of each HLA-II pocket as well as a general view cut through the HLA-II groove
	
Both can be run on any copy of the asymmetric unit by signalling the chains of the complex using --chains. 
The chains argument must be in the order DRA, DRB, peptide; typically (and for these structures, always) --chains ABC or --chains DEF etc.

## Getting Started

To run both analyses on every structure, use:
* analyse_everything.sh

This should generate all structural images, providing the .pdb and .mtz files are downloaded, named accordingly and located within the working directory.

The bash file will demonstrate how to setup arguments, however, use --help on either script to get information about arguments for each.
Both scripts dependent on modules defined in bin/ and well as data info in bin/data/

If you wish to run the analysis on your own structures, please add your peptide name to the dictionary in bin/data/peptide_colours.py and then point towards this colour using the --peptide argument. 

## Dependencies

The following packages are required in your python environment (tested on Python 3.8.13 - Anaconda distribution - linux-64)
* pymol-open-source
 
 The scripts make a number of subprocess.call requests and therefore require the following in your path, all of which are achieved through installing CCP4 or Phenix or both:
 
 Phenix:
* phenix.refine
* phenix.mtz2map
* phenix.pdbtools
 
 CCP4

* refmac5
* pdbcur
* mapmask
* fft
 
 For omit map analysis there are two options. Use pMHC_analyse_v2.py --omit_mode [phenix/refmac] to select. Refmac is faster, Phenix is slower. Phenix mode may be more robust as it will implement simmulated annealing after atom omission.

I intend to neaten up, build upon and generalise this analysis as a pipeline in the future..
-- Bruce