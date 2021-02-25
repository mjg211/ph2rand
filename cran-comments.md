## Test environments
* local Windows 10, R 4.0.4
* local OS X 11.12.1 install, R 4.0.3

## R CMD check results
There were no ERRORs or WARNINGs.

There was 1 NOTE on OS X:

* checking installed package size ...
     installed size is 12.8Mb
     sub-directories of 1Mb or more:
       doc   10.8Mb
       help   1.2Mb

There was 1 NOTE on Windows 10:

* checking compiled code ... NOTE
  Note: information on .o files for i386 is not available
  Note: information on .o files for x64 is not available
  File '~/ph2rand.Rcheck/ph2rand/libs/i386/ph2rand.dll':
    Found 'abort', possibly from 'abort' (C), 'runtime' (Fortran)
    Found 'exit', possibly from 'exit' (C), 'stop' (Fortran)
    Found 'printf', possibly from 'printf' (C)
  File '~/ph2rand.Rcheck/ph2rand/libs/x64/ph2rand.dll':
    Found 'abort', possibly from 'abort' (C), 'runtime' (Fortran)
    Found 'exit', possibly from 'exit' (C), 'stop' (Fortran)
    Found 'printf', possibly from 'printf' (C)

  Compiled code should not call entry points which might terminate R nor
  write to stdout/stderr instead of to the console, nor use Fortran I/O
  nor system RNGs. The detected symbols are linked into the code but
  might come from libraries and not actually be called.

  See 'Writing portable packages' in the 'Writing R Extensions' manual.

## Downstream dependencies
There are currently no downstream dependencies for this package.