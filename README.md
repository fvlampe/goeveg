# goeveg R-package
Functions for Community Data and Ordinations

A collection of functions useful in (vegetation) community analyses and ordinations. The ordination functions work as an addition to the functions from the `vegan`-package. 

Includes:
* Automatic species selection for ordination diagrams based on cover abundances and species fit (`ordiselect` - function)
* Generation of species response curves (`specresponse` - function)
* Scree/stress plots for NMDS (`screeplot_NMDS` - function)
* Rank-abundance curve plotting for single or multiple samples (`racurve` and `racurves`-functions).
* Calculation and sorting of synoptic tables with fidelity and differential species assessment (`syntable` and `synsort` functions)
* Taxa merging for taxa with identical names (`merge_taxa` - function)
* One-step cleaning and transposing of vegetation matrices: (`clean_matrix` and `trans_matrix` - functions)
* Conversion between cover-abundance codes and percentage cover (`cov2per` and `per2cov` - functions)

Furthermore some basic functions are included, such as standard error of the mean `sem`, coefficient of variance `cv` or conversion between degrees and radians `deg2rad`/`rad2deg`.

