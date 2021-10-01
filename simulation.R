##This script is to test ssGSEA (implemented in the GSVA package) in single cell RNAseq analysis.
##See paper https://www.biorxiv.org/content/10.1101/2021.06.29.450404v2
##Created by Siyuan Zheng (zhengs3@uthscsa.edu) and Nighat Noureen (noureen@uthscsa.edu)
##Last modified: 10/1/2021

##why should we be cautious about using ssGSEA for single cell RNAseq analysis?
##The method was first described by Barbie et al. (Nature, 2009). 
##ssGSEA calcuates signature scores by summing up differences in ECDF between genes in and out of a signature. 
##This substantially  differs from the original GSEA algorithm (Subramanian et al. PNAS 2005) that seeks max differences between the two ECDFs.
##The ssGSEA scoring scheme relies on gene ranks. While these ranks are largely continuous in bulk expression data, they can be disrupted in single cell data due to dropouts. Dropout genes have an expression value of zero, creating tie ranks for all dropout genes.  


library(GSVA)
library(gplots)

#----------------------------------------------------------------------------------------------------------#
#Step I: generate an expression matrix, with columns representing samples and rows representing genes
#Gene expression is ranked from high to low, the same as described in the original paper (Barbie et al. 2009)
#The matrix has 10,000 genes, 9 dropout rates ranging from 0 to 0.8 increased by 0.1
#Each column/sample has the exact same expression pattern except differences in dropouts.

N=10000  #total number of genes
M=N:1    ##initial expression matrix, corresponding to zero dropout 
dropout.rate=seq(0.1,0.8,0.1)
for (rate in dropout.rate){
	n_zero=N*rate
	value=1:N
	value[1:n_zero]=0
	value=sort(value,decreasing=T)
	M=cbind(M,value)
}
rownames(M)=paste('g',1:N,sep='')
colnames(M)=paste('drop',c(0,dropout.rate),sep='') 

#----------------------------------------------------------------------------------------------------------#
#Step II: generate random gene sets, and compare their scores at different dropout rates.
##Each gene signature is 100 genes for convenience
##Signatures can mimic different scenarios, such as top genes plus "noises", evenly distributed genes, etc. 
##To get a robust pattern, signatures of the same distribution is generated for n.iteration times.
##note that, once a gslist is generated, jump to run/visual step. 

n.iteration=1000

#senario I: random sampling from the entire gene pool, equivalent to uniform distribution. 
gslist=list()
for (i in 1:n.iteration){
	gs=sample(rownames(M),100)
	gslist[[i]]=gs
}

#senario II, signature comprises mostly top genes, with few noises
gslist=list()
for (i in 1:n.iteration){
	gs1=sample(rownames(M)[1:3000],80)
	gs2=sample(rownames(M)[5001:10000],20)
	gslist[[i]]=c(gs1,gs2)
}

#senario III, signature comprises mostly top genes, with few noises
gslist=list()
for (i in 1:n.iteration){
	gs1=sample(rownames(M)[1:5000],80)
	gs2=sample(rownames(M)[5001:7000],20)
	gslist[[i]]=c(gs1,gs2)
}

##senario IV, signature genes comparises mostly down genes, with few noises
gslist=list()
for (i in 1:n.iteration){
	gs1=sample(rownames(M)[1:3000],20)
	gs2=sample(rownames(M)[4001:8000],80)
	gslist[[i]]=c(gs1,gs2)
}


#----------------------------------------------------------------------------------------------------------#
##Step III, run ssGSEA and visualize outputs.
##For visualization, we use the sample without dropout as reference. Blue--score decrease; Red--score increase. 

result=gsva(M,gset.idx.list=gslist,method='ssgsea',ssgsea.norm=F)
x11(width=8,height=8)
delta=apply(result,2,function(x) x-result[,1])
heatmap.2(delta,Colv=F,Rowv=F,trace='none',density='none',scale='none',col=bluered(20))


#----------------------------------------------------------------------------------------------------------#
#We build a toy example to further demonstrate the impact of dropouts on ssGSEA scores. We use a signature comprising n+1 genes.
#We first randomly select 99 genes from the top 1000 genes of the dummy expression matrix. These 99 genes are fixed. 
##Then, starting from the 1001 genes onwards, we add one gene to the signature each time, starting from the 1001st gene throughout the list. 
##This will allow us to examine how the signature score changes as the additional gene gradually moves from a high expression position to dropouts. 

start.gs=sample(rownames(M)[1:1000],99)
gslist=list(start.gs)
for (i in 1001:10000){
	newgs=c(start.gs, paste('g',i,sep=''))
	gslist=c(gslist,list(newgs))
}
result=gsva(M,gset.idx.list=gslist,method='ssgsea',ssgsea.norm=F)

lowerb=floor(min(result)/50)*50
upperb=ceiling(max(result)/50)*50

plot(result[,'drop0'],axes=F,ylim=c(lowerb,upperb),ylab='ssGSEA score') ##scores based on the original profile (no dropout, models bulk samples)
axis(1,labels=seq(1000,10000,1000),at=seq(0,9000,1000))
axis(2,at=seq(lowerb,upperb,25))
points(result[,'drop0.1'],col='gold') #scores based on 10% dropout
points(result[,'drop0.3'],col='red')  #scores based on 30% dropout
points(result[,'drop0.6'],col='blue') #scores based on 6% dropout
points(result[,'drop0'],col='black')  #scores are identical for the non-dropout zone


