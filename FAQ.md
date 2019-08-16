#### I want to ransform GFF file using sam file, what parameter could I tune?
The -l and -w are important. Simply speaking, if you set them with larger values, the result should be better and it would take more computational time.

#### When I run gean, I got "Illegal instruction" error
Some AVX2 code was used to speed up the ZSDP algorithm implemented GEAN. If you got this error, this suggest you CPU does not support AVX2 technology. <br />
We suggested you to change another computer to run GEAN, since some advanced functions only work with AVX2 supporting. <br />
