ClonalTREE (Clonal reconstruction in Time couRse long-term Evolution Experiment)  

Usage: python3 ClonalTREE.py \<VAF file\> \<algorithm\> \<fail threshold\> \<out prefix\> \<clones (optional)\>  

VAF file:       [String] Input file containing the variant allele frequencies matrix (F).  
algorithm:      [Int] 0 (RP-RT); 1 (EP-ET); 2 (EP-GT); 3 (GP-GT).  
fail threshold: [Float] Minimum value allowed in the C matrix to define failure.  
out prefix:     [String] File path to prefix all output files.  
clones:         [String] A file containing the composition of known clones. (Optional argument)  
  
Example VAF File contents:  
0 1 2 3 4 5 6 7 8 9  
0.94968 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000  
0.70349 0.31292 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000  
0.38256 0.08753 0.38538 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000  
0.21387 0.08375 0.57983 0.04710 0.00000 0.00000 0.00000 0.00000 0.00000  
0.26783 0.19797 0.48729 0.11363 0.26040 0.08971 0.00000 0.00000 0.00000  
0.12287 0.11063 0.33476 0.03762 0.04898 0.39338 0.11046 0.00000 0.00000  
0.13605 0.13088 0.22169 0.06396 0.11725 0.50164 0.24372 0.13126 0.10530  

Example clones file contents:  
1 2  
6 7 8  


