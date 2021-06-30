################################################################################################################################################################################

### JASMINE (Jointly Assessing Signature Mean and INferring Enrichment):- Single Cell RNAseq signature scoring Method using set of marker genes

####   Input Data:- 
####                 1. A single cell RNAseq data matrix, where rows correspond to gene symbols and columns correspond to cells
####                 2. A vector of marker genes(symbols) reflecting a biological process for scoring

####   Output :-
######                Single Cells with corresponding scores for biological process.


###     Authors:-  Nighat Noureen and Siyuan Zheng
###     Email:-   noureen@uthscsa.edu
###               ZhengS3@uthscsa.edu



############# How to Call JASMINE Function ###########

## data = Input data
### genes = List of gene symbols for which scores have to be calculated

Result  =  JASMINE(data,genes) ## calling JASMINE 



#################################################################################################################################################################################
## libraries and functions required

stringsAsFactors=FALSE
library(stringr)
library(GSA)

##################################### STRUCTURE of JASMINE ####################
### Function1:-  Calculating Mean Ranks for signature genes across each cell

RankCalculation <- function(x,genes){
			
            subdata = x[x!=0]                                                                      ### Removing Dropouts from single cell
			DataRanksUpdated=rank(subdata)                                                         ### Calculating ranks of each signature gene per cell
			DataRanksSigGenes = DataRanksUpdated[which(names(DataRanksUpdated) %in% genes)]        ### Shortling rank vector for signature genes 
			CumSum = ifelse(length(DataRanksSigGenes),mean(DataRanksSigGenes,na.rm = TRUE),0 )     ### Calculating Mean of ranks for signature genes
 			FinalRawRank = CumSum/length(subdata)                                                  ### Normalizing Means by total coverage
			return(FinalRawRank)                                                
			}			

#### Function2:- Calculating enrichment of signature genes across each cell 	

ORCalculation <- function(data,genes){
			GE = data[which(rownames(data) %in% genes),]                                          ### Subsetting data for signature genes
 			NGE = data[-which(rownames(data) %in% genes),]                                        ### Subsetting data for non-signature genes
			SigGenesExp = apply(GE,2,function(x) length(x[x!=0]))                                 ### Calculating Number of expressed Signature Genes per cell
			NSigGenesExp = apply(NGE,2,function(x) length(x[x!=0]))                               ### Calculating Number of expressed Non-Signature Genes per cell
			SigGenesNE = nrow(GE) - SigGenesExp                                                   ### Calculating Number of Not expressed Signature Genes per cell
			SigGenesNE = replace(SigGenesNE,SigGenesNE==0,1)									  ### Replacing Zero's with 1
			NSigGenesExp = replace(NSigGenesExp,NSigGenesExp==0,1)                                ### Replacing Zero's with 1
		    NSigGenesNE = nrow(data) - (NSigGenesExp + SigGenesExp)                               ### Calculating Number of Not expressed Non-Signature Genes per cell
			NSigGenesNE = NSigGenesNE - SigGenesNE
			OR = (SigGenesExp * NSigGenesNE) / (SigGenesNE * NSigGenesExp)                         ### Calculating Enrichment (Odds Ratio)
            return(OR)
			}

###  Function 3:- Scalar [0,1] Normalization of Means and Enrichment across set of cells

NormalizationJAS <- function(JAS_Scores)
            {
				JAS_Scores = (JAS_Scores - min(JAS_Scores))/(max(JAS_Scores)- min(JAS_Scores))
				return(JAS_Scores)
			}


### Function 4:- Signature Scoring via JASMINE mergining Means and Enrichment
			
JASMINE <- function(data,genes)
		{
  		    idx = match(genes,rownames(data))                                                 
	        idx = idx[!is.na(idx)]
			if(length(idx)> 1){
			RM = apply(data,2,function(x) RankCalculation(x,genes))                              ### Mean RankCalculation for single cell data matrix
			OR = ORCalculation(data,genes)			                                             ### Signature Enrichment Calculation for single cell data matrix
			RM = NormalizationJAS(RM)                                                            ### Normalizing Mean Ranks
            OR = NormalizationJAS(OR)															 ### Normalizing Enrichment Scores
			JAS_Scores = (RM + OR)/2
			FinalScores = data.frame(names(RM),JAS_Scores)                                       ### JASMINE scores
			colnames(FinalScores)[1]='SampleID'
			return(FinalScores)
			}
		}

