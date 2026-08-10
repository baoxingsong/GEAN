#### I want to transform a GFF file using a SAM file. Which parameters can I tune?

The `-l` and `-w` parameters are the most important. In general, increasing their values should improve the results, but it will also increase the computational time.

#### When I run GEAN, I get an "Illegal instruction" error

GEAN uses AVX2 instructions to accelerate its implementation of the ZSDP algorithm. This error usually indicates that your CPU does not support AVX2.

We recommend running GEAN on a computer with AVX2 support because some advanced features require it.
