# This folder contains scripts develop to support collaboration with Michael Grigg's Lab at NIAID.

## get_allele_freq_from_allele.pl

get_allele_freq_from_allele.pl -i < allele file >
 
Optional parameters (defaults indicated within “[]”):

-c <minimum read coverage of allele [2]

-aa <print only alt alleles T/F [F]

-q <mean read quality [10]

-rb <min fwd/rev read ratio [0.1]

-rt <max fwd/rev read ratio [0.9]>

### Example: 
```
get_allele_freq_from_allele.pl -i Mixedchcpcm.allele -aa T -q 60 -rb 0.3 -rt 0.7 -c 30 > my_output.txt
```
