## Test Environments

* local Windows 10 install, R 4.1
* local CentOS 7 install, R 4.1
* Windows-latest (on GitHub action), R-release
* Ubuntu 20.04 (on GitHub action), R-release and R-devel

## CRAN Test Information

Proper use of this package requires a Julia and DifferentialEquations.jl installation.
This is noted in the installation guide. The Github Actions tests show that on
the major operating systems, if this is installed, then the package will successfully
pass its tests. However, since these softwares are not available on all of the CRAN
computers the tests fail as expected there, and are thus skipped.

Previously, the package was archived on CRAN due to accidentally installing
Julia on the CRAN computers. Given the comments from the Twitter thread
https://twitter.com/ChrisRackauckas/status/1405518810958991361, the community
suggested the fix of using \dontrun instead of \donttest to the offending examples, 
which is the fix that has been implemented here.

## R CMD check results

I would like to keep the CITATION.bib file to conform to the standard format
so that way it can get indexed by non-R tools. The R standard citation is
also included in inst.

## Downstream Dependencies

N/A

## Authors and Copyright

The copyright is held by the SciML organization, which the sole author Chris
Rackauckas is on the steering committee. The other contributions were deemed to
not substantial to merit authorship, being at most ~40 characters total.
