# goeveg R-package
Functions for Community Data and Ordinations

A collection of functions useful in (vegetation) community analyses and ordinations. The ordination functions work as add-on for functions from `vegan`-package. 

Includes:
* Automatic selection of species for ordination diagrams using limits for cover abundances and/or species fit  (`ordiselect()` - function)
* Species response curves (`specresponse()` - function)
* Scree/stress plots for NMDS (`screeplot_NMDS()` - function)
* Rank-abundance curves (`racurve()` and `racurves()` - functions; the latter for (multiple) samples)
* Calculation and sorting of synoptic tables with calculation of fidelity and differential species (`syntable()` and `synsort()` functions)
* Merging of taxa with identical names (`merge_taxa` - function)
* Conversion between cover-abundance codes and percentage coverages (`cov2per` and `per2cov` - functions)


Furthermore some basic functions are included, such as standard error of the mean `sem()`, coefficient of variance `cv()` or conversion between degrees and radians `deg2rad`/`rad2deg`.

