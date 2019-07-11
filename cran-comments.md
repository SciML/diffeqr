## Test Environments

* local Windows 10 install, R 3.4.4
* local CentOS 7 install, R 3.4.4
* Travis-CI Ubuntu 14.04 LTS, R 3.4.4
* Travis-CI MacOSX 16.7.0, R 3.4.4
* AppVeyor Windows Server 2012 R2 x64 (build 9600), R 3.4.4

## CRAN Test Information

Proper use of this package requires a Julia and DifferentialEquations.jl installation.
This is noted in the installation guide. The Travis and AppVeyor tests show that on
the major operating systems, if this is installed, then the package will successfully
pass its tests. However, since these softwares are not available on all of the CRAN
computers the tests fail as expected there, and are thus skipped.

## R CMD check results

I would like to keep the CITATION.bib file to conform to the standard format
so that way it can get indexed by non-R tools. The R standard citation is
also included in inst

## Downstream Dependencies

N/A
