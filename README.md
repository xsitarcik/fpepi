# FPepi: Flower pollination algorithm for detection of epistasis associated with a phenotype
The FPepi tool uses flower pollination algorithm coupled with the tabu search to find potential SNP combinations, that are then evaluated by the G-test. FPepi has two variants, the first uses two objective functions: K2-score and Gini score, while the second uses K2-score and mutual information score. Objectives are optimized separately using two populations.

#### Requirements
The FPepi tool requires Python 3.6 with numpy, pandas and scipy library. 

### Input parameters
Following parameters, that are required to be set inside each script, are as follows:
`gini_population` represents the number of agents in a population of flowers with Gini score as its objective    
`k2_population` represents the number of agents in a population of flowers with K2 score as its objective    
`iteration_size` denotes the number of iterations for which the FPepi algorithm will be running     
`alpha_value` denotes p-value threshold (before Bonferroni correction) which must be passed by the SNP combination to be said as significant    
`searching_path` path to the input file     
`output_file` path to the file, where the results of the epiBAT algorithm will be outputted    

There are also optional parameters, that can be modified:  
`init_prob` represents the initial probability of switch probability in flower pollination algorithm    
`min_dof` the minimum required number of samples in one column of contingency table, in other case, the column is discarded    
`index_beta` beta parameter of Levy flight    
`required_iterations_for_taboo` represents the number of iterations when the best solution in a population is not changed, that is needed to add that best solution to the tabu table    
`zeta_radius` defines the approximity when comparing solutions    
`bee_limit` represents the number of iterations when a solution is not changed, that is then discarded    
`reflection_coef` represents reflection coefficient of simplex method    
`extension_coef` represents extension coefficient of simplex method    
`compression_coef` represents compression coefficient of simplex method    
`shrink_coef` represents shrink coefficient of simplex method    
`elite_crossover` represents elite crossover paramater, i.e. probability that one SNP of SNP combination will be replaced by other SNP randomly selected from an agent of population     

# How to run FPepi
### Input data file format
Input data must be in a usual comma-delimited file containing the case-control genotype data. The first line in a file denotes the SNP IDs, whereas the last column denotes the class, i.e. case or control. The following rows contain the genotype data and the disease state, while the genotype data should have values 0, 1, 2 (i.e. the homozygous major allele, heterozygous allele, and homozygous minor allele), and the disease state should have values 0 and 1 (i.e. control and case). 
Example of the input data containing 5 SNPs and 3 samples is as follows:

X0,X1,X2,X3,X4,X5,Class      
0,0,0,1,0,0      
0,2,0,2,1,0      
1,0,0,1,1,1      
  
After setting the input parameters and preparing your input data, open the command line, go to the directory, where you have downloaded the FPepi_gini.py or the FPepi_mi.py script and call:   
`FPepi_mi.py` or `FPepi_gini.py`

FPepi_mi.py is the first variant using K2 score and mutual information score as objectives. FPepi_gini.py is the second variant using K2 score and gini score as objectives. Results shown that gini score achieves high detection power mainly in the disease models with low heritability.

### Output of the FPepi tool
In a specified file, the FPepi saves results of its 1st stage (i.e. candidate set before evaluation by the G-test) and also final results (i.e. candidate set after evaluation by the G-test) with the corresponding p-values.
