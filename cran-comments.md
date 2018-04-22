## Test Environments

* local Windows 10 install, R 3.4.4
* local CentOS 7 install, R 3.4.4
* Travis-CI Ubuntu 14.04 LTS, R 3.4.4
* Travis-CI MacOSX 16.7.0, R 3.4.4
* AppVeyor Windows Server 2012 R2 x64 (build 9600), R 3.4.4

## R CMD check results

2 Errors: Proper use of this package requires a Julia installation. This is noted in the installation
          guide. The Travis and AppVeyor tests show that on the major operating systems, if this is
          installed, then the package will successfully pass its tests.
2 Note: The appveyor.yml and cran-comments.md (this) file are for testing and registering purposes.
        The other note I'm not sure what it's for.

## Downstream Dependencies

N/A
