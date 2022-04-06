# JASMINE
Single cell RNAseq Signature scoring
JASMINE is an R script which could be used directly by using the functions without installation
####   Input Data for JASMINE:- 
####                 1. A single cell RNAseq data matrix, where rows correspond to gene symbols and columns correspond to cells
####                 2. A vector of marker genes(symbols) reflecting a biological process for scoring

####   Output of JASMINE :-
######                Single Cells with corresponding scores for biological process.


###     Authors:-  Nighat Noureen and Siyuan Zheng
###     Email:-   noureen@uthscsa.edu
###               ZhengS3@uthscsa.edu



############# How to Call JASMINE Function ###########

## data = Input data
### genes = List of gene symbols for which scores have to be calculated
### method = Enter the method for computing the enrichment. Method could be either 'oddsratio' or 'likelihood'

Result  =  JASMINE(data,genes,method =c('oddsratio','likelihood')) ## calling JASMINE 
