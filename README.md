# FPepi: Flower pollination algorithm for detection of epistasis associated with a phenotype
The FPepi tool uses flower pollination algorithm coupled with the tabu search to find potential SNP combinations, that are then evaluated by the G-test. FPepi has two variants, the first uses two objective functions: K2-score and Gini score, while the second uses K2-score and mutual information score. Objectives are optimized separately using two populations.

#### Requirements
The FPepi tool requires Python 3.6 with numpy, pandas and scipy library. 

### Input parameters
There are following parameters, that can be set:
TBA

# How to run FPepi
### Input data file format
Input data must be in a usual comma-delimited file containing the case-control genotype data. The first line in a file denotes the SNP IDs, whereas the last column denotes the class, i.e. case or control. The following rows contain the genotype data and the disease state, while the genotype data should have values 0, 1, 2 (i.e. the homozygous major allele, heterozygous allele, and homozygous minor allele), and the disease state should have values 0 and 1 (i.e. control and case). 
Example of the input data containing 5 SNPs and 3 samples is as follows:

X0,X1,X2,X3,X4,X5,Class      
0,0,0,1,0,0      
0,2,0,2,1,0      
1,0,0,1,1,1      
  
After setting the input parameters and preparing your input data, open the command line, go to the directory, where you have downloaded the FPepi.py script and call:   
`FPepi.py`

### Output of the FPepi tool
In a specified file, the FPepi saves results of its 1st stage (i.e. candidate set before evaluation by the G-test) and also final results (i.e. candidate set after evaluation by the G-test) with the corresponding p-values.
