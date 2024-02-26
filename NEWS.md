# goeveg 0.7.4
* Added option to transform with individual scales per sample in 'cov2per' and 'per2cov'
* Added descriptions and option 'rmchar' to 'trans_matrix()'

# goeveg 0.7.3
* Added function 'clean_matrix()' to remove species without occurrences (frequency = 0) and samples without species from a species matrix in one simple step
* Added function 'trans_matrix()' to transpose a species matrix, while preserving correct species and sample names. 
* Simplified function 'merge_taxa()' to make it work much faster.
    * Added option 'backtransform' to decide whether cover-abundance values should be kept as percentage cover or back-transformed into original cover-abundance values
    * Option 'drop_zero' renamed to 'clean_matrix' and set on FALSE by default
* Fixed an error in 'cov2per' when providing a data frame with only one column as community data
* Fixed an unnecessary warning message in 'syntable' occurring at cover value transformation

# goeveg 0.7.2
* Added cover abundance scale "niwap" from Lower Saxony species survey programmes (Schacherer 2001)
* Added 'x'-value to presence/absence scale
* 'merge_taxa()':
    * returns names of merged taxa only once (not for each relevé)
    * added option 'drop_zero' to decide whether species without occurrences or empty samples should be removed or kept
    * fixed an error when providing individual scales for each sample
* 'syntable()' and 'merge_taxa()' automatically repair imported tables with empty character values ("")

# goeveg 0.7.1
* 'syntable()' can now handle factorial variables for defining clusters, e.g. to summarize relevés according to pre-defined categories and is more flexible regarding the format of the community matrix 
* Terminology was harmonized between different functions

# goeveg 0.7.0
* New functions added:
    * 'merge_taxa()' for merging taxa/species with identical names
    * 'cov2per()' and 'per2cov()' for conversion between cover-abundance codes and percentage cover
* 'dimcheckMDS()' is renamed into 'screeplot_NMDS()' with enhanced description and progress bar

# goeveg 0.6.5
* Fixed wrong species labeling in 'racurve()' when 'freq = TRUE'

# goeveg 0.6.4
* Added functions 'deg2rad()' and 'rad2deg()' for conversion between radians and degrees
* Updated data table 'schedenenv'

# goeveg 0.6.3
* Fixed wrong p-value calculation for GLMs in 'specresponse()'
* Fixed problem with NAs in 'specresponse()' when showing point values

# goeveg 0.6.2
* Added na.action argument to 'specresponse()'

# goeveg 0.6.1
* Explained deviances and p-values are now printed in 'specresponse()'. Full model results are returned in an (invisible) object. 
* Added functionality to select the least abundant (rarest) species in 'ordiselect()'

# goeveg 0.6.0
* (Re-)added functions for calculation and sorting of synoptic tables: 'syntable()' and 'synsort()'
* Comprehensive update for 'ordiselect()'. Now returns exact proportion of selected species. Correction in selection to axis fit limits. Variable fit now only works with factor centroids. 
* Updated help pages

# goeveg 0.5.1
* Fixes in references and value tags

# goeveg 0.5.0
* Removal of functions 'synsort()' and 'syntable()' due to unsolved incompatibilities

# goeveg 0.4.4
* Added lwd argument for 'specresponse()'

# goeveg 0.4.3
* Added xlim & ylim arguments for 'racurve()'
* Added na.rm argument for 'ordiselect()'

# goeveg 0.4.2
* Small fixes, fixed package dependencies
* Spell checking

# goeveg 0.4.1
* Added new functions for calculation and sorting of synoptic tables: 'syntable()' and 'synsort()'

# goeveg 0.3.3
* Fixes in documentation

# goeveg 0.3.2
* Merged 'specresponses()'/'specresponse()' into one single function 'specresponse()'
* Better selection method of polynomial GLMs and GAMs in 'specresponse()'

# goeveg 0.3.1

* Fixed use of external functions ('gam()', 'rdist()')
* Max. of 3 polynomials in automatic GLM selection of 'specresponse()'

# goeveg 0.3.0

* Fixed and renewed function 'specresponse()'/'specresponses()': now works also with NMDS, includes zero values and is based on presence/absence data (logistic regression). Instead of cubic smoothing splines the function now uses GLMs/GAMs.

# goeveg 0.2.0

* Added functionality to use frequencies in 'racurve()'
* Added functionality to label species in 'racurve()'
* New (invisible) output in 'racurve()'
* Package checked and tested on OS X

# goeveg 0.1.6

* Fi