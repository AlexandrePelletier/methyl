
renv::install("yuanzhongshang/MRAID")

out<-"outputs/P01-MRAID_eQTL"
dir.create(out)

source("scripts/utils/new_utils.R")
library(data.table)
library(MRAID)


#Sumamry data with 100 candidate correlated SNPs, n1=30000, n2=30000, 
library(MRAID)
data_path<-"tools/MRAID/example/"

#load the Zscore vector for the exposure
x<-fread(fp(data_path,"Zscore_1.txt"))
Zscore_1<-x$x

#load the Zscore vector for the outcome
y<-fread(fp(data_path,"Zscore_2.txt"))
Zscore_2<-y$x


#load the LD matrix in exposure GWAS data 
zx<-fread(fp(data_path,"Sigma1sin.txt"))
Sigma1sin<-as.matrix(data.frame(zx,row.names = "V1"))

#load the LD matrix in outcome GWAS data 
zy<-fread(fp(data_path,"Sigma2sin.txt"))
Sigma2sin<-as.matrix(data.frame(zy,row.names = "V1"))

#load the sample size
samplen1=30000
samplen2=30000

#run MRAID 
?MRAID
result<-MRAID(Zscore_1, Zscore_2, Sigma1sin, Sigma2sin,
              samplen1, samplen2, 
              Gibbsnumber=1000,burninproportion=0.2,
              pi_beta_shape=0.5,
pi_beta_scale=4.5,pi_c_shape=0.5,pi_c_scale=9.5,pi_1_shape=0.5,pi_1_scale=1.5,pi_0_shape=0.05,pi_0_scale=9.95,option=1)

result

#results

#plot - plot effect

