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

* First CRAN release
