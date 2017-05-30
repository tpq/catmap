## catmap 1.6.1.9000
---------------------
* Add unit tests for core functions
    * `test_catmap.R` unit tests complete and validated
    * `test_forest.R` unit tests complete and validated
* Major revisions to `catmap.funnel`
    * Function now rebuilt using `metafor` package
* Major revisions to forest plots
    * New `makeForest` function builds all forest figures
    * Modified `catmap.forest` to use `makeForest`
    * [`catmap.cumulative`]
    * [`catmap.sense`]

## catmap 1.6.1
---------------------
* Major revisions to global package architecture
    * Convert Rd files to Roxygen format
    * Add DOIs to package DESCRIPTION
    * Add NEWS and README files
* Major revisions to `catmap` function
    * Function now takes `data.frame` or `matrix` input
    * Function now returns output as a list object
    * Removed redundant and superfluous code
    * Removed `print.all` argument entirely
    * Revised documentation
* Major revisions to `catmapdata` data
    * Store `catmapdata` object in RDA format
    * Fix WARNING associated with `catmapdata` manual
    * Revised documentation
