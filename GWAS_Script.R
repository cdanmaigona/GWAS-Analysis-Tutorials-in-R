#R packages for GWAS analysis
###install.packages('SNPassoc')
###install.packages('IntAssoPlot')
###install.packages('GenABEL')
###install.packages('rrBLUP')
###BiocManager::install("GWASTools")
###install.packages('GenABEL')
###install.packages('statgenGWAS')
###library(devtools) install_github("covaruber/sommer")
#estimate additive (A.mat), 
#dominance (D.mat), and 
#epistatic (E.mat) relationship matrices to model covariances among genotypes

###https://link.springer.com/protocol/10.1007/978-1-4939-6682-0_14#copyrightInformation 
###LOAD THE WORKSPACE with the file "GWAS.DATA.RData" downloaded from Springer in R. Click in Session and then click on load workspace.
install.packages("rrBLUP")
install.packages("corrplot")

library(rrBLUP)
library(corrplot)

#Load data or read in data
pheno <- read.csv('phenoat.csv',header=T)

#########################################################################
##########################*****PHENOTYPES*****###########################
#########################################################################
head(pheno)### View phenotypic data.
str(pheno) ### Check variables structure. GID and ENV needs to be factors.


###Checking normality of the data. In theory residuals needs to be checked 
###but in general if data are normal, residuals are normal. 

hist(pheno$Yield, col="black",xlab="Yield",ylab="Frequency",
     border="white",breaks=10,main="Yield Histogram") # Data seems pretty normal with some outliers. 

shapiro.test(pheno$Yield)##Shapiro-Wilk test indicates that normality condition is not met. 



## Let´s remove the outliers and check the normality again. 
boxplot.yield <- boxplot(pheno$Yield,xlab="Boxplot",col='purple',ylab="Yield",ylim=c(4000,9000))

outliers <- boxplot.yield$out; outliers #10 outliers

pheno <- pheno[-which(pheno$Yield%in%outliers),] #removing outliers. 


shapiro.test(pheno$Yield)# After removing outliers, Shapiro-Wilk test indicates Normality.
## When data are nor normal, you can improve normality of the original phenotypic data using
# logit, square root, arcsine, and log transformation methods.  


pheno <- na.omit(pheno)## To remove all posible misiing data.

#########################################################################
##########################*****GENOTYPES*****###########################
#########################################################################
####Other R Packeges for SNP data analysis
#Ape
#Adegenet
#Apparent
#Dartr
#snpReady
#genotypeR
#SNPRelate
#vcfR
#Poppr
#phytools

geno[1:5,1:5] ### View genotypic data.
map[1:5,1:3] ### View Map data.

##########################***FILTERING***###############################

## Filtering conditions will depend on the researcher. In this simple, function
## we remove individuals with missing data, markers with a certain % of missing,
## and heterozygous calls (IF there are a high proportion of heterozygous, it can
## indicate a problem with the marker, because oat is an inbreed specie)

filter.fun <- function(geno,IM,MM,H){
  #Remove individuals with more than a % missing data
  
individual.missing <- apply(geno,1,function(x){
    return(length(which(is.na(x)))/ncol(geno))
  })
  
  
  #length(which(individual.missing>0.40)) #will tell you how many 
  #individulas needs to be removed with 20% missing.
  
#Remove markers with % missing data
  marker.missing <- apply(geno,2,function(x)
  {return(length(which(is.na(x)))/nrow(geno))
    
  })
  
  
  length(which(marker.missing>0.6))
  #Remove markers herteozygous calls more than %. 
  heteroz <- apply(geno,1,function(x){
    return(length(which(x==0))/length(!is.na(x)))
  })
  
  filter1 <- geno[which(individual.missing<IM),which(marker.missing<MM)]
  filter2 <- filter1[,(heteroz<H)]
  return(filter2)
}
geno[1:10,1:10]

geno.filtered <- filter.fun(geno[,1:3629],0.4,0.60,0.02)
geno.filtered[1:5,1:5];dim(geno.filtered)


##########################***IMPUTATION***###############################

### rrBLP program will make imputation. For the simplicity , we impute using 
# the mean but EM algorithm can be also used. 
## rrBLUP also allows to remove markers depending on the Minor allele frequency (MAF),
## in our example we remove those markers with MAF less than 0.05.

library(rrBLUP)
Imputation <- A.mat(geno.filtered,impute.method="EM",return.imputed=T,min.MAF=0.05)

K.mat <- Imputation$A ### KINSHIP matrix
geno.gwas <- Imputation$imputed #NEW geno data.
geno.gwas[1:5,1:5]## view geno
K.mat[1:5,1:5]## view Kinship


################***CHECKING POPULATION STRUCTURE EFFECTS***###############
## Principal components analysis

geno.scale <- scale(geno.gwas,center=T,scale=F) # Data needs to be center.
svdgeno <- svd(geno.scale) 
PCA <- geno.scale%*%svdgeno$v #Principal components
PCA[1:5,1:5]
## Screeplot to visualize the proportion of variance explained by PCA

plot(round((svdgeno$d)^2/sum((svdgeno$d)^2),d=7)[1:10],type="o",main="Screeplot",xlab="PCAs",ylab="% variance")

##Proportion of variance explained by PCA1 and PCA2
PCA1 <- 100*round((svdgeno$d[1])^2/sum((svdgeno$d)^2),d=3); PCA1
PCA2 <- 100*round((svdgeno$d[2])^2/sum((svdgeno$d)^2),d=3); PCA2
PCA3 <- 100*round((svdgeno$d[3])^2/sum((svdgeno$d)^2),d=3); PCA3
PCA4 <- 100*round((svdgeno$d[3])^2/sum((svdgeno$d)^2),d=3); PCA4
### Plotting Principal components.

plot(PCA[,1],PCA[,2],xlab=paste("Pcomp:",PCA1,"%",sep=""),ylab=paste("Pcomp:",PCA2,"%",sep=""),pch=20,cex=0.7)

### Plotting depending on clustering. 
Eucl <- dist(geno.gwas) ###Euclinean distance
fit <- hclust(Eucl,method="ward.D2")### Ward criterion makes clusters with same size.
groups2 <- cutree(fit,k=3) ### Selecting two clusters.
table(groups2)# Number of individuals per cluster.

plot(PCA[,1],PCA[,2],xlab=paste("Pcomp:",PCA1,"%",sep=""),ylab=paste("Pcomp:",PCA2,"%",sep=""),pch=20,cex=0.7,col=groups2)
legend("bottomright",c("Subpop1: 244","Subpop2: 84"),pch=20,col=(c("black","red")),lty=0,bty="n",cex=1)

################***MATCHING PHENOTYPE AND GENOTYPE***###############

pheno=pheno[pheno$GID%in%rownames(geno.gwas),]
pheno$GID<-factor(as.character(pheno$GID), levels=rownames(geno.gwas)) #to assure same levels on both files
pheno <- pheno[order(pheno$GID),]
##Creating file for GWAS function from rrBLUP package
X<-model.matrix(~-1+ENV, data=pheno)
pheno.gwas <- data.frame(GID=pheno$GID,X,Yield=pheno$Yield)
head(pheno.gwas)

geno.gwas <- geno.gwas[rownames(geno.gwas)%in%pheno.gwas$GID,]
pheno.gwas <- pheno.gwas[pheno.gwas$GID%in%rownames(geno.gwas),]
geno.gwas <- geno.gwas[rownames(geno.gwas)%in%rownames(K.mat),]
K.mat <- K.mat[rownames(K.mat)%in%rownames(geno.gwas),colnames(K.mat)%in%rownames(geno.gwas)]
pheno.gwas <- pheno.gwas[pheno.gwas$GID%in%rownames(K.mat),]

################***MATCHING GENOTYPE AND MAP***###############
geno.gwas<-geno.gwas[,match(map$Markers,colnames(geno.gwas))]
head(map)
geno.gwas <- geno.gwas[,colnames(geno.gwas)%in%map$Markers]
map <- map[map$Markers%in%colnames(geno.gwas),]
geno.gwas2<- data.frame(mark=colnames(geno.gwas),chr=map$chrom,loc=map$loc,t(geno.gwas))
dim(geno.gwas2)
colnames(geno.gwas2)[4:ncol(geno.gwas2)] <- rownames(geno.gwas)

head(pheno.gwas)
geno.gwas2[1:6,1:6]
K.mat[1:5,1:5]

#####################################***ANALYSIS***################################################
gwasresults<-GWAS(pheno.gwas,geno.gwas2, fixed=colnames(pheno.gwas)[2:5], K=NULL, plot=T,n.PC=0)
gwasresults2<-GWAS(pheno.gwas,geno.gwas2, fixed=colnames(pheno.gwas)[2:5], K=NULL, plot=T,n.PC=6)
gwasresults3<-GWAS(pheno.gwas,geno.gwas2, fixed=colnames(pheno.gwas)[2:5], K=K.mat, plot=T,n.PC=0)
gwasresults4<-GWAS(pheno.gwas,geno.gwas2, fixed=colnames(pheno.gwas)[2:5], K=K.mat, plot=T,n.PC = 6)

#The option plot=T will produce manhattan plots and q-q plots. With the aim to provide
#another option in R, another graphs in R have been provided. 
###################################################################################
####################*** QQ-MANHATTAN and CORRELATION PLOTS*****####################
###################################################################################
#Let´s see the structure
str(gwasresults)
str(gwasresults)
#First 3 columns are just the information from markers and map.
#Fouth and next columns are the results form GWAS. Those values are already
#the  -log10 pvalues, so no more transformation needs to be done to plot them. 

###################################################################################
#################################*** QQ PLOT*****##################################
###################################################################################


  N <- length(gwasresults$Yield)
  expected.logvalues <- sort( -log10( c(1:N) * (1/N) ) )
  observed.logvalues <- sort( gwasresults$Yield)
  
  plot(expected.logvalues , observed.logvalues, main="Naïve model(K=NULL,n.PC=0)", 
       xlab="expected -log pvalue ", 
       ylab="observed -log p-values",col.main="blue",col="coral1",pch=20)
  abline(0,1,lwd=3,col="black")
  
  
  N1 <- length(gwasresults2$Yield)
  expected.logvalues1 <- sort( -log10( c(1:N1) * (1/N1) ) )
  observed.logvalues1 <- sort( gwasresults2$Yield)
  
  plot(expected.logvalues1 , observed.logvalues1, main="Q model (K=NULL,n.PC=6)", 
       xlab="expected -log pvalue ", 
       ylab="observed -log p-values",col.main="blue",col="coral1",pch=20)
  abline(0,1,lwd=2,col="black")
  
  
  N2 <- length(gwasresults3$Yield)
  expected.logvalues2 <- sort( -log10( c(1:N2) * (1/N2) ) )
  observed.logvalues2 <- sort( gwasresults3$Yield)
  
  plot(expected.logvalues2 , observed.logvalues2, main="K model (K=Kmat,n.PC=0)", 
       xlab="expected -log pvalue ", 
       ylab="observed -log p-values",col.main="blue",col="coral1",pch=20)
  abline(0,1,lwd=2,col="black")
  
  N3 <- length(gwasresults4$Yield)
  expected.logvalues3 <- sort( -log10( c(1:N3) * (1/N3) ) )
  observed.logvalues3 <- sort( gwasresults4$Yield)
  
  plot(expected.logvalues3 , observed.logvalues3, main="Q+K model (K.mat,n.PC=6)", 
       xlab="expected -log pvalue ", 
       ylab="observed -log p-values",col.main="blue",col="coral1",pch=20)
  abline(0,1,lwd=2,col="black")
  

  
  ###################################################################################
  #################################*** MANHATTAN PLOT*****###########################
  ###################################################################################
  #False Discovery Rate Function
  
  FDR<-function(pvals, FDR){
    pvalss<-sort(pvals, decreasing=F)
    m=length(pvalss)
    cutoffs<-((1:m)/m)*FDR
    logicvec<-pvalss<=cutoffs
    postrue<-which(logicvec)
    print(postrue)
    k<-max(c(postrue,0))
    cutoff<-(((0:m)/m)*FDR)[k+1]
    return(cutoff)
  }
  
  alpha_bonferroni=-log10(0.05/length(gwasresults$Yield)) ###This is Bonferroni correcton
  alpha_FDR_Yield <- -log10(FDR(10^(-gwasresults$Yield),0.05))## This is FDR cut off
  
  
  #################################*** MANHATTAN PLOT*****###########################
  
  plot(gwasresults$Yield,col=gwasresults$chr,ylab="-log10.pvalue",
       main="Naïve model (K=NULL,n.PC=0)",xaxt="n",xlab="Position",ylim=c(0,14))
  #axis(1,at=c(1:length(unique(gwasresults$chr))),labels=unique(gwasresults$chr))
  axis(1,at=c(0,440,880,1320,1760))
  abline(a=NULL,b=NULL,h=alpha_bonferroni,col="blue",lwd=2)
  abline(a=NULL,b=NULL,h=alpha_FDR_Yield,col="red",lwd=2,lty=2)
  legend(1,13.5, c("Bonferroni","FDR") , 
         lty=1, col=c('red', 'blue'), bty='n', cex=1,lwd=2)
  
  plot(gwasresults2$Yield,col=gwasresults2$chr,ylim=c(0,14),ylab="-log10.pvalue",
       main="Q model (K=NULL,n.PC=6)",xaxt="n",xlab="Position")
  axis(1,at=c(0,440,880,1320,1760))
  abline(a=NULL,b=NULL,h=alpha_bonferroni,col="blue",lwd=2)
  abline(a=NULL,b=NULL,h=alpha_FDR_Yield,col="red",lwd=2,lty=2)
  legend(1.5,13.5, c("Bonferroni","FDR") , 
         lty=1, col=c('red', 'blue'), bty='n', cex=1,lwd=2)
  
  plot(gwasresults3$Yield,col=gwasresults3$chr,ylim=c(0,14),ylab="-log10.pvalue",
       main="K model (K=K.mat,n.PC=0)",xaxt="n",xlab="Position")
  axis(1,at=c(0,440,880,1320,1760))
  abline(a=NULL,b=NULL,h=alpha_bonferroni,col="blue",lwd=2)
  abline(a=NULL,b=NULL,h=alpha_FDR_Yield,col="red",lwd=2,lty=2)
  legend(1.5,13.5, c("Bonferroni","FDR") , 
         lty=1, col=c('red', 'blue'), bty='n', cex=1,lwd=2)
  
  plot(gwasresults4$Yield,col=gwasresults4$chr,ylim=c(0,14),ylab="-log10.pvalue",
       main="Q+K model (K=K.mat,n.PC=6)",xaxt="n",xlab="Position")
  axis(1,at=c(0,440,880,1320,1760))
  abline(a=NULL,b=NULL,h=alpha_bonferroni,col="blue",lwd=2)
  abline(a=NULL,b=NULL,h=alpha_FDR_Yield,col="red",lwd=2,lty=2)##FDR gives inf for Yield
  legend(1,13.5, c("Bonferroni","FDR") , 
         lty=1, col=c('red', 'blue'), bty='n', cex=1,lwd=2)

  
  ############################*** WHICH ARE HITS?*****###########################
  
  which(gwasresults$Yield>alpha_bonferroni)
  which(gwasresults$Yield>alpha_FDR_Yield)
  which(gwasresults2$Yield>alpha_bonferroni)
  which(gwasresults2$Yield>alpha_FDR_Yield)
  which(gwasresults3$Yield>alpha_bonferroni)
  which(gwasresults3$Yield>alpha_FDR_Yield)
  which(gwasresults4$Yield>alpha_bonferroni)
  which(gwasresults4$Yield>alpha_FDR_Yield)
  
  markers.gwasresults4.bonf<- geno.gwas[,c(53,56,57,1054,1427)]#gwasresults3 and 4 have same hits.
  markers.gwasresults2.bonf <- geno.gwas[,c(53,56,57,1054,1427,1428)]
  
  ###################################################################################
  ##############################*** CORRELATION PLOT*****############################
  ###################################################################################
  
  corr_sign <- cor(markers.gwasresults4.bonf,use="complete.obs") 
  
  corrplot(corr_sign, order="hclust", method="pie", tl.pos="lt", type="upper",        
           tl.col="black", tl.cex=0.8, tl.srt=55, 
           sig.level=0.90,cl.length=21,insig = "blank")     
  mtext("Correlation Significant hits",outer=TRUE,line=1) 
 
#############################Bootstrapping & Phylogeny###############################################
#############################Tutorial and script available for free on github########################
#############################https://github.com/elinck/molecologist##################################

#References
#Isidro-S�nchez J., Akdemir D., Montilla-Basc�n G. (2017) Genome-Wide Association Analysis Using R. 
#In: Gasparis S. (eds) Oat. Methods in Molecular Biology, vol 1536. Humana Press, New York, NY. https://doi.org/10.1007/978-1-4939-6682-0_14








