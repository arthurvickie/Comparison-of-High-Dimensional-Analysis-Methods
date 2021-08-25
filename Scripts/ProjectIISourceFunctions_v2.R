############################################
## Source functions for Project II       ###
############################################
###################################   
DetermineAssocRSNPs = function(gene, LowLD, percentageAssoc){
  #determine which SNPs to actually set as assoc. based on gene
  #check gene
  if(gene == "NAT2"){
    if(LowLD == TRUE){
      #then looking at SNPs in low LD
      if(percentageAssoc == 5){
        #then 1 SNP is associated
        assocSNPs = c(1)
      } else if(percentageAssoc == 15){
        #then 2 SNPs are associated
        assocSNPs = c(1,2)
      } else {
        #then 5 SNPs are associated
        assocSNPs = c(1,2,9,7,10)
      } 
    } else {
      #looking at SNPs in high LD 
      if(percentageAssoc == 5){
        #then 1 SNP is associated
        assocSNPs = c(13)
      } else if(percentageAssoc == 15){
        #then 2 SNPs are associated
        assocSNPs = c(13,14)
      } else {
        #then 5 SNPs are associated
        assocSNPs = c(13,14,4,5,3)
      } 
    }
  } else if(gene == "CHI3L2"){
    if(LowLD == TRUE){
      #then looking at SNPs in low LD
      if(percentageAssoc == 5){
        #then 2 SNPs are associated
        assocSNPs = c(12,2)
      } else if(percentageAssoc == 15){
        #then 5 SNPs are associated
        assocSNPs = c(12,2,3,9,8)
      } else {
        #then 8 SNPs are associated
        assocSNPs = c(12,2,3,9,8,24,28,19)
      } 
    } else {
      #looking at SNPs in high LD 
      #then looking at SNPs in low LD
      if(percentageAssoc == 5){
        #then 2 SNPs are associated
        assocSNPs = c(17,29)
      } else if(percentageAssoc == 15){
        #then 5 SNPs are associated
        assocSNPs = c(5,6,23,18,13)
      } else {
        #then 8 SNPs are associated
        assocSNPs = c(5,6,23,18,13,15,17,29)
      }
    }
  } else { 
    #otherwise we have ASAH1
    if(LowLD == TRUE){
      #then looking at SNPs in low LD
      if(percentageAssoc == 5){
        #then 2 SNPs are associated
        assocSNPs = c(31,9)
      } else if(percentageAssoc == 15){
        #then 6 SNPs are associated
        assocSNPs = c(31,9,4,8,10,5)
      } else if(percentageAssoc == 25){
        #then 10 SNPs are associated
        assocSNPs = c(31,9,4,8,10,5,3,2,7,1)
      }
    } else {
      #looking at SNPs in high LD 
      #then looking at SNPs in low LD
      if(percentageAssoc == 5){
        #then 2 SNPs are associated
        assocSNPs = c(25, 32)
      } else if(percentageAssoc == 15){
        #then 6 SNPs are associated
        assocSNPs = c(37,29,36,28,30,33)
      } else if(percentageAssoc == 25){
        #then 10 SNPs are associated
        assocSNPs = c(37,29,36,28,30,33,34,38,25,32)
      } 
    }
  }
  return(assocSNPs)
}  

##Step 1) Generate paired genotype data, get this data into correct format
#sample haplotypes from 1000 Genomes
sampleHaplotypes = function(chr = c(), i = c(), numPairs = c(), RAND = c(), gene="", path=paste0("/path/to/data/"), outfile = paste0("outfile.vcf")){
  #function to sample Haplotypes from 1K genomes
  #writes vcf file to output, doesn't directly return anything
  
  #make sure number of pairs and random number are integer types
  numPairs = as.integer(numPairs)
  RAND = as.integer(RAND)
  
  #setwd to location of the data
  setwd(path)
  
  #read in marker information from Haploview
  markerInfo=read.table(paste0("AllMarkers_",gene,"_MAF05.txt"), header = F) 
  
  #remove SNPs with MAF less than 0.05 and with percent genotyped over 75 (less than 25% missing)
  #and with HW p values less than 0.001
  snpSubset=which((markerInfo[,7]>=75)&(markerInfo[,10]>=0.05)&(markerInfo[,6]>0.001))  
  
  #subset the data for the VCF file
  markerSubset = markerInfo[snpSubset,]
  
  #read in the haplotype frequency file
  hapFreqFile=read.table(paste0("HapFreqAllMarkers_",gene,"_MAF05_spaced.txt"), header = F)
  
  ###Minor allele frequencies
  allMinorAlleles=substring(markerInfo[,11],3,3) #this lists the minor alleles
  #these are the minor alleles as A,G,T,C values for VCF files
  minorAlleleSubset.letters = allMinorAlleles[snpSubset]
  allMinorAlleles[allMinorAlleles=="A"]=1 #recodes the MAs as numbers, so that they match the haplotype freq file
  allMinorAlleles[allMinorAlleles=="C"]=2
  allMinorAlleles[allMinorAlleles=="G"]=3
  allMinorAlleles[allMinorAlleles=="T"]=4
  minorAlleleSubset=as.numeric(allMinorAlleles)[snpSubset] #this lists the minor alleles as numbers for all SNPs that pass filtering
  
  #this lists the major alleles
  allMajorAlleles = substring(markerInfo[,11],1,1) 
  #this lists the major alleles as numbers for all SNPs that pass filtering
  majorAlleleSubset.letters = allMajorAlleles[snpSubset] 
  
  #vectorizes the haplotype frequencies and combines them
  haplotypeFreq.vec=as.vector(t(hapFreqFile[,1:(ncol(hapFreqFile)-3)]))
  
  #checks if haplotypeFreq.vec values are equal to MA, basically like, yes minor allele = 1, no = 0
  #dimension is number of haplotypes possible x number of SNPs in haplotype
  ma.haplotypeFreq.vec = matrix(as.numeric(haplotypeFreq.vec==rep(minorAlleleSubset,nrow(hapFreqFile))),nrow=nrow(hapFreqFile),ncol=length(minorAlleleSubset),byrow=T)
  
  #divides haplotype frequencies by total sum (which is slightly less than 1) to get accurate freq values
  correctedHapFreqs = hapFreqFile[,ncol(hapFreqFile)-1] /sum(hapFreqFile[,ncol(hapFreqFile)-1])
  
  #not sure what KK stands for...
  numDonors = numRecip = numPairs
  KK = numDonors + numRecip
  
  #should put in a random seed here to make sure sample differs each time
  set.seed(RAND)
  
  #pulling indices of the haplotype frequencies based on prob of haplotype, pulls KK of them
  idHap1 = sample(1:length(correctedHapFreqs),KK,replace = T,prob = correctedHapFreqs) 
  #pulling indices of the haplotype frequencies based on prob of haplotype, pulls KK of them
  idHap2 = sample(1:length(correctedHapFreqs),KK,replace = T,prob = correctedHapFreqs) 
  #this is adding the 2 haplotypes to make full genotype (I think I use this?)
  #file is length of haplotype (number of SNPs) x KK in dimension
  sampledGenotype = ma.haplotypeFreq.vec[idHap1, ] + ma.haplotypeFreq.vec[idHap2, ] 
  
  ############################################################################################
  ##make these genotypes into VCF files
  #Then VCF -> plink
  
  #create vcf template that will get filled in by the needed information, but contains all the regular info
  vcfTemplate = matrix(NA, nrow=length(minorAlleleSubset), ncol=(9 + KK))
  #nrow is number of SNPs, ncols is  9 + number of subjects generated
  
  #first col should be chrom number
  vcfTemplate[,1] = chr
  
  #second column is position, need to fill in later
  
  #third column is id, which is just a . for msprime vcfs, 6 (quality) and 8 (info) also have .s
  vcfTemplate[,3] = "."
  vcfTemplate[,6] = "."
  vcfTemplate[,8] = "."
  
  #4th column is reference allele, and 5th is alternative allele, msprime vcf has A and T
  vcfTemplate[,4] = majorAlleleSubset.letters
  vcfTemplate[,5] = minorAlleleSubset.letters
  
  #7th column is filter, msprime vcf has PASS
  vcfTemplate[,7] = "PASS"
  
  #9th column is format, msprime vcf has GT for genotype
  vcfTemplate[,9] = "GT"
  
  #add the vcf template to the final list
  finalVCFsVar = vcfTemplate
  
  ###################################################
  #now work on adding the positions to the vcfs
  finalVCFsVar[,2] = markerSubset$V3
  ###################################################
  #adding genotypes to the final vcfs
  
  #make new lists to put output in
  hapsSubsetVarMerged = matrix(paste(t(ma.haplotypeFreq.vec[idHap1,]), t(ma.haplotypeFreq.vec[idHap2,]), sep="|"), nrow=length(minorAlleleSubset), ncol=(KK))
  
  #add to the final vcfs file
  finalVCFsVar[,10:(KK+9)] = hapsSubsetVarMerged
  
  ##########################################################
  #now need to export VCF file, sans comments, as separate files
  write.table(finalVCFsVar, file=outfile, sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
}

#pull R genotype information
obtainRGenotypes = function(chr = c(), numSamples = c(), simNum = c(), gene = "", path=paste0("/path/to/data/")){
  #function of obtain Recipient genotypes for D/R transplant pairs
  #output is a matrix of R genotypes for all individuals (N x m)
  #load needed packages
  suppressMessages(require(ARTP2,lib.loc='/home/vlynn/R/library'))
  suppressMessages(require(dplyr))
  
  #make sure number of samples is integer value
  numSamples = as.integer(numSamples) 
  
  #setwd to location of plink data
  setwd(path)
  
  #make lists of names of plink data
  #for HapGen generated data
  bedFile = paste0(gene,"_",numSamples,"Pairs_SimNum",simNum,"_Subset.bed")
  bimFile = paste0(gene,"_",numSamples,"Pairs_SimNum",simNum,"_Subset.bim")
  famFile = paste0(gene,"_",numSamples,"Pairs_SimNum",simNum,"_Subset.fam")
  
  #read in plink binary files
  plinkFile = read.bed(bed = bedFile, bim = bimFile, fam = famFile)
  
  #matched pairs are even and odd columns of each file
  #define variables for number of columns and rows of the df
  nrowPlink = nrow(plinkFile) #this is number of generated subjects (2*SS)
  ncolPlink = ncol(plinkFile) #this is number of SNPs
  
  recipientGenotypes = plinkFile[seq(2,nrowPlink,2),]
  RGenosMat = matrix(unlist(recipientGenotypes), ncol = ncolPlink, byrow = F)
  return(RGenosMat)
}
#pull D genotype information
obtainDGenotypes = function(chr = c(), numSamples = c(), simNum = c(), gene = "", path=paste0("/path/to/data/")){
  #function of obtain Donor genotypes for D/R transplant pairs
  #output is a matrix of D genotypes for all individuals (m x N)
  #load needed packages
  suppressMessages(require(ARTP2,lib.loc='/home/vlynn/R/library'))
  suppressMessages(require(dplyr))
  
  #make sure number of samples is integer value
  numSamples = as.integer(numSamples) 
  
  #setwd to location of plink data
  setwd(path)
  
  #make lists of names of plink data
  #for HapGen generated data
  bedFile = paste0(gene,"_",numSamples,"Pairs_SimNum",simNum,"_Subset.bed")
  bimFile = paste0(gene,"_",numSamples,"Pairs_SimNum",simNum,"_Subset.bim")
  famFile = paste0(gene,"_",numSamples,"Pairs_SimNum",simNum,"_Subset.fam")
    
  #read in plink binary files
  plinkFile = read.bed(bed = bedFile, bim = bimFile, fam = famFile)
  
  #matched pairs are even and odd columns of each file
  #define variables for number of columns and rows of the df
  nrowPlink = nrow(plinkFile) #this is number of generated subjects (2*SS)
  ncolPlink = ncol(plinkFile) #this is number of SNPs
  
  #Ds are odd numbered rows, Rs are even numbered rows
  donorGenotypes = plinkFile[seq(1,nrowPlink,2),]
  DGenosMat = matrix(unlist(donorGenotypes), ncol = ncolPlink, byrow = F)
  return(DGenosMat)
  
}

###################################       
##Step 2) Calculate individual scores for each pair,
##then calc score for gene region
#use R and D genotypes to calc IBS mismatch score
calcIBSMismatch = function(RGenosMat = matrix(), DGenosMat = matrix()){
  #function to calculate the IBS mismatch score for all D/R pairs
  #returns a matrix of IBS mismatch scores  (m x N)
  
  #number of SNPs is the number of columns in D Genotypes matrix
  ncolPlink = ncol(DGenosMat)

  #calculate the difference between the two subjects
  diffsPlink = abs(DGenosMat - RGenosMat)
  
  #######################
  ## IBS Score
  ####################### 
  #if diff = 0, score = 0
  #if diff = 1, score is unchanged
  #if diff = 2, score = 2
  IBSMismatch = diffsPlink
  
  #save IBS scores
  RIBSScoresMat = matrix(unlist(IBSMismatch), ncol = ncolPlink, byrow = F)
  return(RIBSScoresMat)
}
#use R and D genotypes to calc Incompatibility score
calcIncompatibilityScore = function(RGenosMat = matrix(), DGenosMat = matrix()){
  #function to calculate the incompatibility score for all D/R pairs
  #returns a matrix of incompatibility scores  (m x N)
  
  #number of SNPs is the number of columns in D Genotypes matrix
  ncolPlink = ncol(DGenosMat)
  
  #calculate the difference between the two subjects
  diffsPlink = abs(DGenosMat - RGenosMat)
  
  #######################
  ## Incomp Score
  #######################
  #initialize a list of empty dfs with same number of columns as original
  incomp = diffsPlink
  for(ii in 1:ncol(diffsPlink)){
    incomp[diffsPlink[,ii] == 0,ii] = 0
  }
  for(ii in 1:ncol(diffsPlink)){
    incomp[diffsPlink[,ii] != 0,ii] = 1
  }
  
  #save incomp scores
  RIncompScoresMat = matrix(unlist(incomp), ncol = ncolPlink, byrow = F)
  return(RIncompScoresMat)
}
#use R and D genotypes to calc AMS score
calcAMS = function(RGenosMat = matrix(), DGenosMat = matrix()){
  #function to calculate the AMS score for all D/R pairs 
  #returns a matrix of AMS scores (m x N)
  
  #number of SNPs is the number of columns in D Genotypes matrix
  ncolPlink = ncol(DGenosMat)
  
  #calculate the difference between the two subjects
  diffsPlink = abs(DGenosMat - RGenosMat)
  
  #######################
  ## AMS
  #######################
  #mismatch if D has allele not in R
  #sum across both alleles in genotype
  #Score is either 0, 1, or 2
  alloMismatch = matrix(0, nrow = nrow(diffsPlink), ncol = ncol(diffsPlink)) #make default value 0
  alloMismatch[(DGenosMat == 0) & (RGenosMat == 2)] = 2 #Donor AA, Recip aa
  alloMismatch[(DGenosMat == 2) & (RGenosMat == 0)] = 2 #Donor aa, Recip AA
  alloMismatch[(DGenosMat == 1) & (RGenosMat == 2)] = 1 #Donor Aa, Recip aa
  alloMismatch[(DGenosMat == 1) & (RGenosMat == 0)] = 1 #Donor Aa, Recip AA
  alloMismatch[is.na(DGenosMat) | is.na(RGenosMat)] = NA #make sure NAs are preserved
  #match row and column names from the original data set
  #row names should be R Ids
  rownames(alloMismatch) = rownames(DGenosMat)
  colnames(alloMismatch) = colnames(DGenosMat)
  
  #save AMS scores
  RAlloMismatchMat = matrix(unlist(alloMismatch), ncol = ncolPlink, byrow = F)
  return(RAlloMismatchMat)
}
#use R and D genotypes to calc binary mismatch score
calcBinaryMM = function(RGenosMat = matrix(), DGenosMat = matrix()){
  #function to calculate the binary MM score for all D/R pairs
  #returns a matrix of binary MM scores  (m x N) 
  
  #number of SNPs is the number of columns in D Genotypes matrix
  ncolPlink = ncol(DGenosMat)
  
  #calculate the difference between the two subjects
  diffsPlink = abs(DGenosMat - RGenosMat)
  
  #######################
  ## Binary  Mismatch
  #######################
  #mismatch if D has allele not in R
  #Score is either 0 or 1
  binMismatch = matrix(0, nrow = nrow(diffsPlink), ncol = ncol(diffsPlink)) #make default value 0
  binMismatch[(DGenosMat == 1) & (RGenosMat == 2)] = 1 #Donor Aa, Recip aa
  binMismatch[(DGenosMat == 0) & (RGenosMat == 2)] = 1 #Donor AA, Recip aa
  binMismatch[(DGenosMat == 2) & (RGenosMat == 0)] = 1 #Donor aa, Recip AA
  binMismatch[(DGenosMat == 1) & (RGenosMat == 0)] = 1 #Donor Aa, Recip AA
  binMismatch[is.na(DGenosMat) | is.na(RGenosMat)] = NA #make sure NAs are preserved
  #match row and column names from the original data set
  #row names should be R Ids
  rownames(binMismatch) = rownames(DGenosMat)
  colnames(binMismatch) = colnames(DGenosMat)
  
  #save Binary mismatch scores
  RBinMismatchMat = matrix(unlist(binMismatch), ncol = ncolPlink, byrow = F)
  return(RBinMismatchMat)
}

#use single SNP kernel functions to calc gene-based score
calcGeneScore = function(SingleSNPKernel = matrix()){
  #function to calculate the gene based score based on the single SNP
  #kernels
  #Returns 1 x N vector, with single gene score for each individual
  geneScore = rowSums(SingleSNPKernel)
  geneScore.mat = as.matrix(geneScore)
  return(geneScore.mat)
}

#use single SNP kernel functions to calc gene-based score, using percent of SNPs
calcGeneScorePercentOfSNPs = function(SingleSNPKernel = matrix(), gene =  "", percentageAssoc = 100, LowLD = TRUE){
  #function to calculate the gene based score based on the single SNP
  #Only uses a percentage of SNPs in the gene region, assuming not all SNPs are associated with outcome
  #kernels and optional weights
  #Returns 1 x N vector, with single gene score for each individual
  #check gene
  if(gene == "NAT2"){
    if(LowLD == TRUE){
      #then looking at SNPs in low LD
      if(percentageAssoc == 5){
        #then 1 SNP is associated
        SingleSNPKernel = as.matrix(SingleSNPKernel[,1],ncol=1)
      } else if(percentageAssoc == 15){
		#then 2 SNPs are associated
		SingleSNPKernel = SingleSNPKernel[,c(1,2)]
	  } else if(percentageAssoc == 25){
        #then 5 SNPs are associated
        SingleSNPKernel = SingleSNPKernel[,c(1,2,9,7,10)]
      } else if(percentageAssoc == 50){
        #then 7 SNPs are associated
        SingleSNPKernel = SingleSNPKernel[,c(1,2,9,7,10,6,8)]
      } else if(percentageAssoc == 75){
        #then 11 SNPs are associated
        SingleSNPKernel = SingleSNPKernel[,c(1,2,9,7,10,6,8,11,12,3,5)]
      } else {
        #otherwise use all SNPs
        SingleSNPKernel = SingleSNPKernel
      }
    } else {
      #looking at SNPs in high LD 
      if(percentageAssoc == 5){
        #then 1 SNP is associated
        SingleSNPKernel = as.matrix(SingleSNPKernel[,13],ncol=1)
      } else if(percentageAssoc == 15){
		#then 2 SNPs are associated
		SingleSNPKernel = SingleSNPKernel[,c(13,14)]
	  } else if(percentageAssoc == 25){
        #then 5 SNPs are associated
        SingleSNPKernel = SingleSNPKernel[,c(13,14,4,5,3)]
      } else if(percentageAssoc == 50){
        #then 7 SNPs are associated
        SingleSNPKernel = SingleSNPKernel[,c(11,12,3,5,4,14,13)]
      } else if(percentageAssoc == 75){
        #then 11 SNPs are associated
        SingleSNPKernel = SingleSNPKernel[,c(7,10,6,8,11,12,3,5,4,14,13)]
      } else {
        #otherwise use all SNPs
        SingleSNPKernel = SingleSNPKernel
      }
    }
  } else if(gene == "CHI3L2"){
    if(LowLD == TRUE){
      #then looking at SNPs in low LD
      if(percentageAssoc == 5){
        #then 2 SNPs are associated
        SingleSNPKernel = SingleSNPKernel[,c(12,2)]
      } else if(percentageAssoc == 15){
		#then 5 SNPs are associated
		SingleSNPKernel = SingleSNPKernel[,c(12,2,3,9,8)]
	  } else if(percentageAssoc == 25){
        #then 8 SNPs are associated
        SingleSNPKernel = SingleSNPKernel[,c(12,2,3,9,8,24,28,19)]
      } else if(percentageAssoc == 50){
        #then 17 SNPs are associated
        SingleSNPKernel = SingleSNPKernel[,c(12,2,3,9,8,24,28,19,21,22,20,11,1,10,16,4,31)]
      } else if(percentageAssoc == 75){
        #then 25 SNPs are associated
        SingleSNPKernel = SingleSNPKernel[,c(12,2,3,9,8,24,28,19,21,22,20,11,1,10,16,4,31,26,25,32,14,27,33,30,7)]
      } else {
        #otherwise use all SNPs
        SingleSNPKernel = SingleSNPKernel
      }
    } else {
      #looking at SNPs in high LD 
      if(percentageAssoc == 5){
        #then 2 SNPs are associated
        SingleSNPKernel = SingleSNPKernel[,c(17,29)]
      } else if(percentageAssoc == 15){
		#then 5 SNPs are associated
		SingleSNPKernel = SingleSNPKernel[,c(5,6,23,18,13)]
	  } else if(percentageAssoc == 25){
        #then 8 SNPs are associated
        SingleSNPKernel = SingleSNPKernel[,c(5,6,23,18,13,15,17,29)]
      } else if(percentageAssoc == 50){
        #then 17 SNPs are associated
        SingleSNPKernel = SingleSNPKernel[,c(31,26,25,32,14,27,33,30,7,5,6,23,18,13,15,17,29)]
      } else if(percentageAssoc == 75){
        #then 25 SNPs are associated
        SingleSNPKernel = SingleSNPKernel[,c(21,22,20,11,1,10,16,4,31,26,25,32,14,27,33,30,7,5,6,23,18,13,15,17,29)]
      } else {
        #otherwise use all SNPs
        SingleSNPKernel = SingleSNPKernel
      }
    }
  } else { 
    #otherwise we have ASAH1
    if(LowLD == TRUE){
      #then looking at SNPs in low LD
      if(percentageAssoc == 5){
        #then 2 SNPs are associated
        SingleSNPKernel = SingleSNPKernel[,c(31, 9)]
      } else if(percentageAssoc == 15){
		#then 6 SNPs are associated
		if(dim(SingleSNPKernel)[2] < 40){
			SingleSNPKernel = SingleSNPKernel[,c(31,9,4,8,10,5)]
		} else {
			SingleSNPKernel = SingleSNPKernel[,c(31,9,4,8,10,40)]
		}
	  } else if(percentageAssoc == 25){
        #then 10 SNPs are associated
		if(dim(SingleSNPKernel)[2] < 40){
			SingleSNPKernel = SingleSNPKernel[,c(31,9,4,8,10,5,3,2,7,1)]
		} else {
			SingleSNPKernel = SingleSNPKernel[,c(31,9,4,8,10,40,5,3,2,7)]
		}
      } else if(percentageAssoc == 50){
        #then 20 SNPs are associated
        SingleSNPKernel = SingleSNPKernel[,c(31,9,4,8,10,40,5,3,2,7,1,35,21,22,11,6,15,12,13,14)]
      } else if(percentageAssoc == 75){
        #then 30 SNPs are associated
        SingleSNPKernel = SingleSNPKernel[,c(31,9,4,8,10,40,5,3,2,7,1,35,21,22,11,6,15,12,13,14,18,39,16,19,20,24,17,26,23,27)]
      } else {
        #otherwise use all SNPs
        SingleSNPKernel = SingleSNPKernel
      }
    } else {
      #looking at SNPs in high LD 
      if(percentageAssoc == 5){
        #then 2 SNPs are associated
        SingleSNPKernel = SingleSNPKernel[,c(25, 32)]
      } else if(percentageAssoc == 15){
		#then 6 SNPs are associated
		SingleSNPKernel = SingleSNPKernel[,c(30,33,34,38,25,32)]
	  } else if(percentageAssoc == 25){
        #then 10 SNPs are associated
        SingleSNPKernel = SingleSNPKernel[,c(27, 37, 29, 36, 28, 30, 33, 34, 38, 25, 32)]
      } else if(percentageAssoc == 50){
        #then 20 SNPs are associated
        SingleSNPKernel = SingleSNPKernel[,c(18,39,16,19,20,24,17,26,23,27,37,29,36,28,30,33,34,38,25,32)]
      } else if(percentageAssoc == 75){
        #then 30 SNPs are associated
        SingleSNPKernel = SingleSNPKernel[,c(1,35,21,22,11,6,15,12,13,14,18,39,16,19,20,24,17,26,23,27,37,29,36,28,30,33,34,38,25,32)]
      } else {
        #otherwise use all SNPs
        SingleSNPKernel = SingleSNPKernel
      }
    }
  }
  
  geneScore = rowSums(SingleSNPKernel)
  geneScore.mat = as.matrix(geneScore)
  return(geneScore.mat)
}

###################################       
##Step 3) Generate covariates
#generate covariate data
GenCovData = function(SampleSize, BinaryValues, ContinuousValues){
  #####################################
  # Simulate covariates
  # returns a matrix of size N x K (N = sample size)
  # where K is number of covariates generated
  # allows input of number of binary and
  # continuous covariate values
  #####################################
  # Binary covariates are drawn from a random
  # binimial distribution with p=0.5
  # Continuous covariate sare drawn from a standard
  # Normal distribution
  
  K = BinaryValues + ContinuousValues
  ss = SampleSize
  
  #define W matrix
  W = matrix(NA,nrow = ss, ncol = K)
  
  for(ii in 1:BinaryValues){
    W[,ii] = rbinom(ss,1,0.5)
  }
  
  for(ii in (BinaryValues+1):K){
    W[,ii] = rnorm(ss, mean = 0, sd = 1)
  }
  
  return(W)
}

###################################       
##Step 4) Generate true phenotype data using covariates only (TIE)
## or using covariates, R genotype, and score (power)
#generate phenotype data - null model, with covariates
GenNullPhenos = function(SampleSize, includeCov = FALSE, YCat = TRUE, YPrev,  Covariates){
  # generates null phenotypes to be used as "true" Y values
  # returns N x 1 vector of phenotypes
  # SampleSize (N) is number of D/R pairs
  # includeCov is FALSE if no covariates are being used, TRUE if they are being used
  # YCat = TRUE if outcome is categorical, otherwise Y is continuous
  # YPrev is the probability that Y=1, only needed if YCat=TRUE
  # Covariates is N x K matrix of covariate values (if they are being used)
  
  ss = SampleSize
  #only define W and eff sizes if we have covariates
  if(includeCov == TRUE){
    W = Covariates
    K = ncol(W)
    #define all effect sizes as 0.5, since this is what most people do in 
    #similar papers for TIE simulations
    eff_sizes = matrix(0.5, nrow = K, ncol = 1)
  }
  if(YCat == TRUE){
    prev = YPrev
  }
  
  #if not including covariates
  if(includeCov == FALSE){
    #if Y is catergorical, and no covariates
    # then Y is pulled from a binom(prev) for each individual
    #logit (pi) = alpha_0
    if(YCat == TRUE){
      nullPhenos = rbinom(ss, 1, prev)
    } 
    #otherwise Y is continuous, and no covs
    #then Y is pulled from a standard normal for each individual
    # basically, Y = epsilon (random error)
    else {
      nullPhenos = rnorm(ss, 0, 1)
    }
  }
  # otherwise we are including covariates
  else{
    #if Y is categorical, and covariates
    #then logit Pr(Y=1) = a0 + a*Covariates
    if(YCat == TRUE){
      #alpha 1 
      alphas = W %*% eff_sizes
      #set alpha0 value to start
      alpha_0 = rep(1,ss)
      #calc alpha_0 + alphas
      lin_predict = alpha_0 + alphas
      #calculate p values for binomials, want colMeans of this to be around prevalence
      prob_case = exp(lin_predict)/(1 + exp(lin_predict))
      #change alpha_0 such that prob of case is close to prevalence
      while((mean(prob_case) > prev) == TRUE) {
        #set alpha_0 value such that prob of case is approx prevalence
        alpha_0 = alpha_0 - 0.01
        #calc alpha_0 + alphas
        lin_predict = alpha_0 + alphas
        #calculate p values for binomials, want colMeans of this to be around prevalence
        prob_case = exp(lin_predict)/(1 + exp(lin_predict))
        #use binomial dist to get phenos
        nullPhenos = rbinom(ss,1,prob_case)
      } 
    }
    #otherwise Y is continuous 
    else {
        #alphas 
        alphas = W %*% eff_sizes
        #epsilon
        epsilon = rnorm(ss,0,1)
        #linear model, added in normally dist error
        lin_predict = alphas + epsilon
        #add to null phenos
        nullPhenos = lin_predict
      }
  }
  nullPhenos = as.matrix(nullPhenos, nrow = ss, ncol = 1)
  return(nullPhenos)
}

#generate phenotype data - alt model, with covariates
GenAltPhenos = function(SampleSize, includeCov = FALSE, YCat = TRUE, YPrev, Covariates, RGenoData = matrix(), ScoreData = matrix(), Betas = matrix(), Gamma = c()){
  # generates "true" phenotypes for power analysis using the alternative
  # hypothesis model
  # returns N x 1 vector of phenotypes
  
  # SampleSize: number of D/R pairs
  # includeCov: TRUE or FALSE for if covariates are included in the modeling
  # YCat: TRUE of FALSE for whether Y is categorical or continuous
  # YPrev: if Y is categorical, the p(y=1)
  # Covariates is K x N matrix of covariate values
  # RGenoData is N x m matrix of R genotype values
  # ScoreData is N x 1 vector of gene-based score values
  # Betas: cov effect sizes for R geno, could all be 0 (m x 1)
  # Gamma: cov effect size for score, could all be 0 (1 x 1)
  
  ss = SampleSize
  #only define W and eff sizes if we have covariates
  if(includeCov == TRUE){
    W = Covariates
    K = ncol(W)
    #define all effect sizes as 0.5, since this is what most people do in 
    #similar papers for TIE simulations
    eff_sizes = matrix(0.5, nrow = K, ncol = 1)
  }
  if(YCat == TRUE){
    prev = YPrev
  }
  
  X = RGenoData
  B = as.matrix(Betas, ncol=1)
  Z = ScoreData
  G = matrix(Gamma, nrow = 1, ncol = 1)
  
  #if not including covariates
  if(includeCov == FALSE){
    #if Y is catergorical, and no covariates
    # then Y is calculated using only XBeta and ZGamma values
    #logit (pi) = alpha_0 + XBeta + ZGamma
    if(YCat == TRUE){
      #beta terms
      betas = X %*% B
      #set alpha0 value to start
      alpha_0 = rep(1,ss)
      #gamma term
      gamma = ScoreData %*% Gamma
      #calc alpha_0 + Betas*SNPGeno + Gamma*Score
      lin_predict = alpha_0 + betas + gamma
      #calculate p values for binomials, want colMeans of this to be around prevalence
      prob_case = exp(lin_predict)/(1 + exp(lin_predict))
      #change beta_0 such that prob of case is close to prevalence
      while((mean(prob_case) > prev) == TRUE) {
        #set beta0 value such that prob of case is approx prevalence
        alpha_0 = alpha_0 - 0.01
        #calc alpha_0 + Betas*SNPGeno + Gamma*Score
        lin_predict = alpha_0 + betas + gamma
        #calculate p values for binomials, want colMeans of this to be around prevalence
        prob_case = exp(lin_predict)/(1 + exp(lin_predict))
        #use binomial dist to get phenos
        altPhenos = rbinom(ss,1,prob_case)
      } 
    } else {
	#otherwise Y is continuous, and no covs
    #then Y is pulled from a standard normal for each individual
    # basically, Y = epsilon (random error)
      #beta terms 
      betas = X %*% B
      #gamma term
      gamma = ScoreData %*% Gamma
      #epsilon
      epsilon = rnorm(ss,0,1)
      #linear model, added in normally dist error
      lin_predict = betas + gamma + epsilon
      #add to null phenos
      altPhenos = lin_predict
    }
  } else { 
  # otherwise we are including covariates
    #if Y is categorical, and covariates
    #then logit Pr(Y=1) = a0 + a*Covariates + XBeta + ZGamma
    if(YCat == TRUE){
      #alpha terms
      alphas = W %*% eff_sizes
      #beta terms
      betas = X %*% B
      #set alpha0 value to start
      alpha_0 = rep(1,ss)
      #gamma term
      gamma = ScoreData %*% G
      #calc alpha_0 + Betas*SNPGeno + Gamma*Score
      lin_predict = alpha_0 + alphas + betas + gamma
      #calculate p values for binomials, want colMeans of this to be around prevalence
      prob_case = exp(lin_predict)/(1 + exp(lin_predict))
      #change beta_0 such that prob of case is close to prevalence
      while((mean(prob_case) > prev) == TRUE) {
        #set beta0 value such that prob of case is approx prevalence
        alpha_0 = alpha_0 - 0.01
        #calc alpha_0 + Betas*SNPGeno + Gamma*Score
        lin_predict = alpha_0 + alphas + betas + gamma
        #calculate p values for binomials, want colMeans of this to be around prevalence
        prob_case = exp(lin_predict)/(1 + exp(lin_predict))
        #use binomial dist to get phenos
        altPhenos = rbinom(ss,1,prob_case)
      } 
    } else {
	  #otherwise Y is continuous 
      #alpha terms
      alphas = W %*% eff_sizes
      #beta terms 
      betas = X %*% B
      #gamma term
      gamma = ScoreData %*% G
      #epsilon
      epsilon = rnorm(ss,0,1)
      #linear model, added in normally dist error
      lin_predict = alphas + betas + gamma + epsilon
      #add to null phenos
      altPhenos = lin_predict
    }
  }
  altPhenos = as.matrix(altPhenos, nrow = ss, ncol = 1)
  return(altPhenos)
}

################################################
## Putting some functions together
CalcNullPhenotypeData = function(chr, numPairs, simNum, YPrev, gene, path, covs){
  #pull recipient and donor genotypes
  RGenos = obtainRGenotypes(chr = chr, numSamples = numPairs, simNum = simNum, gene = gene, path = path)
  DGenos = obtainDGenotypes(chr = chr, numSamples = numPairs, simNum = simNum, gene = gene, path = path)
  
  #calculate single snp scores
  IBS.snp = calcIBSMismatch(RGenosMat = RGenos, DGenosMat = DGenos)
  Incomp.snp = calcIncompatibilityScore(RGenosMat = RGenos, DGenosMat = DGenos)
  AMS.snp = calcAMS(RGenosMat = RGenos, DGenosMat = DGenos)
  BinMM.snp = calcBinaryMM(RGenosMat = RGenos, DGenosMat = DGenos)
  
  #calculate gene based scores
  IBS.gene = calcGeneScore(SingleSNPKernel = IBS.snp)
  Incomp.gene = calcGeneScore(SingleSNPKernel = Incomp.snp)
  AMS.gene = calcGeneScore(SingleSNPKernel = AMS.snp)
  BinMM.gene = calcGeneScore(SingleSNPKernel = BinMM.snp)
  
  if(covs == FALSE){
    #generate phenotypes, both continuous and binary
    CatPhenos = GenNullPhenos(SampleSize = numPairs, includeCov = FALSE, YCat = TRUE, YPrev = YPrev)
    ContPhenos = GenNullPhenos(SampleSize = numPairs, includeCov = FALSE, YCat = FALSE)
    
    allData = list(RGenos, IBS.gene, Incomp.gene, AMS.gene, BinMM.gene, CatPhenos, ContPhenos)
  } else {
    #generate covariates
    #for now, a single binary and a single continous covariate
    CovData = GenCovData(SampleSize = numPairs, BinaryValues = 1, ContinuousValues = 1)
    
    #generate phenotypes, both continuous and binary
    CatPhenos = GenNullPhenos(SampleSize = numPairs, includeCov = TRUE, YCat = TRUE, YPrev = YPrev,  Covariates = CovData)
    ContPhenos = GenNullPhenos(SampleSize = numPairs, includeCov = TRUE, YCat = FALSE,  Covariates = CovData)
    
    allData = list(RGenos, IBS.gene, Incomp.gene, AMS.gene, BinMM.gene, CovData, CatPhenos, ContPhenos)
  }
  return(allData)
}

CalcAltPhenotypeData_Scores = function(chr, numPairs, simNum, YPrev, gene, path, covs = FALSE, percentageAssoc, LowLD, TrueScore){
  if(ORSize == "small"){
    Gamma = c(0.14)
  } else if(ORSize == "medium"){
    Gamma = c(0.41)
  } else {
    Gamma = c(0.69)
  }
  
  #pull recipient and donor genotypes
  RGenos = obtainRGenotypes(chr = chr, numSamples = numPairs, simNum = simNum, gene = gene, path = path)
  DGenos = obtainDGenotypes(chr = chr, numSamples = numPairs, simNum = simNum, gene = gene, path = path)
  
  #calculate single snp scores
  IBS.snp = calcIBSMismatch(RGenosMat = RGenos, DGenosMat = DGenos)
  Incomp.snp = calcIncompatibilityScore(RGenosMat = RGenos, DGenosMat = DGenos)
  AMS.snp = calcAMS(RGenosMat = RGenos, DGenosMat = DGenos)
  BinMM.snp = calcBinaryMM(RGenosMat = RGenos, DGenosMat = DGenos)
  
  #calculate gene based scores
  IBS.gene = calcGeneScore(SingleSNPKernel = IBS.snp)
  Incomp.gene = calcGeneScore(SingleSNPKernel = Incomp.snp)
  AMS.gene = calcGeneScore(SingleSNPKernel = AMS.snp)
  BinMM.gene = calcGeneScore(SingleSNPKernel = BinMM.snp)
  
  #need to use TrueScore to pull gene based scores matrix for generating phenotypes
  if(TrueScore == "IBS"){
    PhenoScore = calcGeneScorePercentOfSNPs(SingleSNPKernel = IBS.snp, gene = gene, percentageAssoc = percentageAssoc, LowLD = LowLD)
  } else if(TrueScore == "Incomp"){
    PhenoScore =   calcGeneScorePercentOfSNPs(SingleSNPKernel = Incomp.snp, gene =  gene, percentageAssoc = percentageAssoc, LowLD = LowLD)
  } else if(TrueScore == "AMS"){
    PhenoScore = calcGeneScorePercentOfSNPs(SingleSNPKernel = AMS.snp, gene =  gene, percentageAssoc = percentageAssoc, LowLD = LowLD)
  } else {
    PhenoScore = calcGeneScorePercentOfSNPs(SingleSNPKernel = BinMM.snp, gene =  gene, percentageAssoc = percentageAssoc, LowLD = LowLD)
  }
  
  #generate covariates
  #for now, a single binary and a single continous covariate
  CovData = GenCovData(SampleSize = numPairs, BinaryValues = 1, ContinuousValues = 1)
  
  #need to define null Betas for phenotype generation
  nSNP = ncol(RGenos) #this should be the number of SNPs
  Betas = rep(0,nSNP) #generate null beta values
  Betas = as.matrix(Betas, ncol = 1)
  
  if(covs == FALSE){
    #if no covariates
    #generate phenotypes, both continuous and binary
    #Based on single true score
    CatPhenos = GenAltPhenos(SampleSize = numPairs, includeCov = FALSE, YCat = TRUE, YPrev = YPrev,  RGenoData = RGenos, ScoreData = PhenoScore, Betas = Betas, Gamma = Gamma)
    ContPhenos = GenAltPhenos(SampleSize = numPairs, includeCov = FALSE, YCat = FALSE,  RGenoData = RGenos, ScoreData = PhenoScore, Betas = Betas, Gamma = Gamma)
    
    AltScorePhenos = list(RGenos, IBS.gene, Incomp.gene, AMS.gene, BinMM.gene, CatPhenos, ContPhenos)
  } else {
    #generate phenotypes, both continuous and binary
    #Based on single true score
    CatPhenos = GenAltPhenos(SampleSize = numPairs, includeCov = TRUE, YCat = TRUE, YPrev = YPrev,  Covariates = CovData, RGenoData = RGenos, ScoreData = PhenoScore, Betas = Betas, Gamma = Gamma)
    ContPhenos = GenAltPhenos(SampleSize = numPairs, includeCov = TRUE, YCat = FALSE,  Covariates = CovData, RGenoData = RGenos, ScoreData = PhenoScore, Betas = Betas, Gamma = Gamma)
    
    AltScorePhenos = list(RGenos, IBS.gene, Incomp.gene, AMS.gene, BinMM.gene, CovData, CatPhenos, ContPhenos)
  }
  return(AltScorePhenos)
}

CalcAltPhenotypeData_RSNPs = function(chr, numPairs, simNum, YPrev, gene, path, covs = FALSE, ORSize, LowLD, percentAssoc, TrueScore){
  #define effect based on OR size
  if(ORSize == "small"){
    effect = 0.14
  } else if(ORSize == "medium"){
    effect = 0.41
  } else {
    effect = 0.69
  }
  
  assocSNPs = DetermineAssocRSNPs(gene = gene, LowLD = LowLD, percentageAssoc = percentageAssoc)
  
  #pull recipient and donor genotypes
  RGenos = obtainRGenotypes(chr = chr, numSamples = numPairs, simNum = simNum, gene = gene, path = path)
  DGenos = obtainDGenotypes(chr = chr, numSamples = numPairs, simNum = simNum, gene = gene, path = path)
  
  #calculate single snp scores
  IBS.snp = calcIBSMismatch(RGenosMat = RGenos, DGenosMat = DGenos)
  Incomp.snp = calcIncompatibilityScore(RGenosMat = RGenos, DGenosMat = DGenos)
  AMS.snp = calcAMS(RGenosMat = RGenos, DGenosMat = DGenos)
  BinMM.snp = calcBinaryMM(RGenosMat = RGenos, DGenosMat = DGenos)
  
  #calculate gene based scores
  IBS.gene = calcGeneScore(SingleSNPKernel = IBS.snp)
  Incomp.gene = calcGeneScore(SingleSNPKernel = Incomp.snp)
  AMS.gene = calcGeneScore(SingleSNPKernel = AMS.snp)
  BinMM.gene = calcGeneScore(SingleSNPKernel = BinMM.snp)
  
  #need to use TrueScore to pull gene based scores matrix for generating phenotypes
  if(TrueScore == "IBS"){
    PhenoScore = IBS.gene
  } else if(TrueScore == "Incomp"){
    PhenoScore = Incomp.gene
  } else if(TrueScore == "AMS"){
    PhenoScore = AMS.gene
  } else {
    PhenoScore = BinMM.gene
  }
  
  #generate covariates
  #for now, a single binary and a single continous covariate
  CovData = GenCovData(SampleSize = numPairs, BinaryValues = 1, ContinuousValues = 1)
  
  #need to define null Betas for phenotype generation
  nSNP = ncol(RGenos) #this should be the number of SNPs
  nullBetas = rep(0,nSNP) #generate null beta values
  Betas = nullBetas
  
  #set assoc Betas
  #all betas have same effect for now
  for(jj in assocSNPs){
    Betas[jj] = effect
    Betas = as.matrix(Betas, ncol = 1)
  }
  
  nullGamma = c(0)
  
  if(covs == FALSE){
    
    #generate phenotypes, both continuous and binary
    #Based on single true score
    CatPhenos = GenAltPhenos(SampleSize = numPairs, includeCov = FALSE, YCat = TRUE, YPrev = YPrev,  RGenoData = RGenos, ScoreData = PhenoScore, Betas = Betas, Gamma = nullGamma)
    ContPhenos = GenAltPhenos(SampleSize = numPairs, includeCov = FALSE, YCat = FALSE,  RGenoData = RGenos, ScoreData = PhenoScore, Betas = Betas, Gamma = nullGamma)
    
    AltPhenos_RSNPs = list(RGenos, IBS.gene, Incomp.gene, AMS.gene, BinMM.gene, CatPhenos, ContPhenos)
  } else {
    
    #generate phenotypes, both continuous and binary
    #Based on single true score
    CatPhenos = GenAltPhenos(SampleSize = numPairs, includeCov = TRUE, YCat = TRUE, YPrev = YPrev,  Covariates = CovData, RGenoData = RGenos, ScoreData = PhenoScore, Betas = Betas, Gamma = nullGamma)
    ContPhenos = GenAltPhenos(SampleSize = numPairs, includeCov = TRUE, YCat = FALSE,  Covariates = CovData, RGenoData = RGenos, ScoreData = PhenoScore, Betas = Betas, Gamma = nullGamma)
    
    AltPhenos_RSNPs = list(RGenos, IBS.gene, Incomp.gene, AMS.gene, BinMM.gene, CovData, CatPhenos, ContPhenos)
  }
  return(AltPhenos_RSNPs)
}

CalcAltPhenotypeData_ScoreAndRSNPs = function(chr, numPairs, simNum, YPrev, gene, path, covs = FALSE, ORSizeGene, ORSizeSNP, LowLDGene, LowLDSNP, percentAssocGene, percentAssocSNP, TrueScore){
  if(ORSizeGene == "small"){
    Gamma = c(0.14)
  } else if(ORSizeGene == "medium"){
    Gamma = c(0.41)
  } else {
    Gamma = c(0.69)
  }
  
  #define effect based on OR size
  if(ORSizeSNP == "small"){
    effect = 0.14
  } else if(ORSizeSNP == "medium"){
    effect = 0.41
  } else {
    effect = 0.69
  }
  
  assocSNPs = DetermineAssocRSNPs(gene = gene, LowLD = LowLDSNP, percentageAssoc = percentAssocSNP)
  
  #pull recipient and donor genotypes
  RGenos = obtainRGenotypes(chr = chr, numSamples = numPairs, simNum = simNum, gene = gene, path = path)
  DGenos = obtainDGenotypes(chr = chr, numSamples = numPairs, simNum = simNum, gene = gene, path = path)
  
  #calculate single snp scores
  IBS.snp = calcIBSMismatch(RGenosMat = RGenos, DGenosMat = DGenos)
  Incomp.snp = calcIncompatibilityScore(RGenosMat = RGenos, DGenosMat = DGenos)
  AMS.snp = calcAMS(RGenosMat = RGenos, DGenosMat = DGenos)
  BinMM.snp = calcBinaryMM(RGenosMat = RGenos, DGenosMat = DGenos)
  
  #calculate gene based scores
  IBS.gene = calcGeneScore(SingleSNPKernel = IBS.snp)
  Incomp.gene = calcGeneScore(SingleSNPKernel = Incomp.snp)
  AMS.gene = calcGeneScore(SingleSNPKernel = AMS.snp)
  BinMM.gene = calcGeneScore(SingleSNPKernel = BinMM.snp)
  
  #need to use TrueScore to pull gene based scores matrix for generating phenotypes
  if(TrueScore == "IBS"){
    PhenoScore = IBS.gene
  } else if(TrueScore == "Incomp"){
    PhenoScore = Incomp.gene
  } else if(TrueScore == "AMS"){
    PhenoScore = AMS.gene
  } else {
    PhenoScore = BinMM.gene
  }
  
  #generate covariates
  #for now, a single binary and a single continous covariate
  CovData = GenCovData(SampleSize = numPairs, BinaryValues = 1, ContinuousValues = 1)
  
  #need to define null Betas for phenotype generation
  nSNP = ncol(RGenos) #this should be the number of SNPs
  nullBetas = rep(0,nSNP) #generate null beta values
  Betas = nullBetas
  
  #set assoc Betas
  #all betas have same effect for now
  for(jj in assocSNPs){
    Betas[jj] = effect
    Betas = as.matrix(Betas, ncol = 1)
  }
  
  nullGamma = c(0)
  
  if(covs == FALSE){
    
    #generate phenotypes, both continuous and binary
    #Based on single true score
    CatPhenos = GenAltPhenos(SampleSize = numPairs, includeCov = FALSE, YCat = TRUE, YPrev = YPrev,  RGenoData = RGenos, ScoreData = PhenoScore, Betas = Betas, Gamma = Gamma)
    ContPhenos = GenAltPhenos(SampleSize = numPairs, includeCov = FALSE, YCat = FALSE,  RGenoData = RGenos, ScoreData = PhenoScore, Betas = Betas, Gamma = Gamma)
    
    AltPhenos_RSNPs = list(RGenos, IBS.gene, Incomp.gene, AMS.gene, BinMM.gene, CatPhenos, ContPhenos)
  } else {
    
    #generate phenotypes, both continuous and binary
    #Based on single true score
    CatPhenos = GenAltPhenos(SampleSize = numPairs, includeCov = TRUE, YCat = TRUE, YPrev = YPrev,  Covariates = CovData, RGenoData = RGenos, ScoreData = PhenoScore, Betas = Betas, Gamma = Gamma)
    ContPhenos = GenAltPhenos(SampleSize = numPairs, includeCov = TRUE, YCat = FALSE,  Covariates = CovData, RGenoData = RGenos, ScoreData = PhenoScore, Betas = Betas, Gamma = Gamma)
    
    AltPhenos_RSNPs = list(RGenos, IBS.gene, Incomp.gene, AMS.gene, BinMM.gene, CovData, CatPhenos, ContPhenos)
  }
  return(AltPhenos_RSNPs)
}

CalcAltPhenotypeData_ScoreAndRSNPs_DiffGenes = function(chrSNPs, chrMatch, numPairs, simNum, YPrev, geneSNPs, geneMatch, path1, path2, covs = FALSE, ORSizeGene, ORSizeSNP, LowLDGene, LowLDSNP, percentAssocGene, percentAssocSNP, TrueScore){
  if(ORSizeGene == "small"){
    Gamma = c(0.14)
  } else if(ORSizeGene == "medium"){
    Gamma = c(0.41)
  } else {
    Gamma = c(0.69)
  }
  
  #define effect based on OR size
  if(ORSizeSNP == "small"){
    effect = 0.14
  } else if(ORSizeSNP == "medium"){
    effect = 0.41
  } else {
    effect = 0.69
  }
  
  assocSNPs = DetermineAssocRSNPs(gene = gene, LowLD = LowLDSNP, percentageAssoc = percentAssocSNP)
  
  #pull recipient and donor genotypes
  RGenosSNPs = obtainRGenotypes(chr = chrSNPs, numSamples = numPairs, simNum = simNum, gene = geneSNPs, path = path1)

  RGenosMatch = obtainRGenotypes(chr = chrMatch, numSamples = numPairs, simNum = simNum, gene = geneMatch, path = path2)
  DGenosMatch = obtainDGenotypes(chr = chrMatch, numSamples = numPairs, simNum = simNum, gene = geneMatch, path = path2)
  
  #calculate single snp scores
  IBS.snp = calcIBSMismatch(RGenosMat = RGenosMatch, DGenosMat = DGenosMatch)
  Incomp.snp = calcIncompatibilityScore(RGenosMat = RGenosMatch, DGenosMat = DGenosMatch)
  AMS.snp = calcAMS(RGenosMat = RGenosMatch, DGenosMat = DGenosMatch)
  BinMM.snp = calcBinaryMM(RGenosMat = RGenosMatch, DGenosMat = DGenosMatch)
  
  #calculate gene based scores
  IBS.gene = calcGeneScore(SingleSNPKernel = IBS.snp)
  Incomp.gene = calcGeneScore(SingleSNPKernel = Incomp.snp)
  AMS.gene = calcGeneScore(SingleSNPKernel = AMS.snp)
  BinMM.gene = calcGeneScore(SingleSNPKernel = BinMM.snp)
  
  #need to use TrueScore to pull gene based scores matrix for generating phenotypes
  if(TrueScore == "IBS"){
    PhenoScore = IBS.gene
  } else if(TrueScore == "Incomp"){
    PhenoScore = Incomp.gene
  } else if(TrueScore == "AMS"){
    PhenoScore = AMS.gene
  } else {
    PhenoScore = BinMM.gene
  }
  
  #generate covariates
  #for now, a single binary and a single continous covariate
  CovData = GenCovData(SampleSize = numPairs, BinaryValues = 1, ContinuousValues = 1)
  
  #need to define null Betas for phenotype generation
  nSNP = ncol(RGenosSNPs) #this should be the number of SNPs
  nullBetas = rep(0,nSNP) #generate null beta values
  Betas = nullBetas
  
  #set assoc Betas
  #all betas have same effect for now
  for(jj in assocSNPs){
    Betas[jj] = effect
    Betas = as.matrix(Betas, ncol = 1)
  }
  
  nullGamma = c(0)
  
  if(covs == FALSE){
    
    #generate phenotypes, both continuous and binary
    #Based on single true score
    CatPhenos = GenAltPhenos(SampleSize = numPairs, includeCov = FALSE, YCat = TRUE, YPrev = YPrev,  RGenoData = RGenosSNPs, ScoreData = PhenoScore, Betas = Betas, Gamma = Gamma)
    ContPhenos = GenAltPhenos(SampleSize = numPairs, includeCov = FALSE, YCat = FALSE,  RGenoData = RGenosSNPs, ScoreData = PhenoScore, Betas = Betas, Gamma = Gamma)
    
    AltPhenos_RSNPs = list(RGenosSNPs, IBS.gene, Incomp.gene, AMS.gene, BinMM.gene, CatPhenos, ContPhenos)
  } else {
    
    #generate phenotypes, both continuous and binary
    #Based on single true score
    CatPhenos = GenAltPhenos(SampleSize = numPairs, includeCov = TRUE, YCat = TRUE, YPrev = YPrev,  Covariates = CovData, RGenoData = RGenosSNPs, ScoreData = PhenoScore, Betas = Betas, Gamma = Gamma)
    ContPhenos = GenAltPhenos(SampleSize = numPairs, includeCov = TRUE, YCat = FALSE,  Covariates = CovData, RGenoData = RGenosSNPs, ScoreData = PhenoScore, Betas = Betas, Gamma = Gamma)
    
    AltPhenos_RSNPs = list(RGenosSNPs, IBS.gene, Incomp.gene, AMS.gene, BinMM.gene, CovData, CatPhenos, ContPhenos)
  }
  return(AltPhenos_RSNPs)
}


RunGLM.LRT = function(allData, FitScore, Cat = TRUE, covs = TRUE){
  #pull elements from allData list
  RGenos = allData[[1]]
  IBS.gene = allData[[2]]
  Incomp.gene = allData[[3]]
  AMS.gene = allData[[4]]
  BinMM.gene = allData[[5]]
  
  if(FitScore == "IBS"){
    Score.gene = IBS.gene
  } else if(FitScore == "Incomp"){
    Score.gene = Incomp.gene
  } else if(FitScore == "AMS"){
    Score.gene = AMS.gene
  } else{
    Score.gene = BinMM.gene
  }
  
  if(covs == FALSE){
    #if no covariates:
    CatPhenos = allData[[6]]
    ContPhenos = allData[[7]]
    
    if(Cat == TRUE){
      #fit the null and alternative models
      fitNull = glm(CatPhenos~RGenos, family = binomial)
      fitAlt = glm(CatPhenos~RGenos+Score.gene, family = binomial, epsilon = 1e-6)
    } else {
      #continuous outcome fit
      fitNull = glm(ContPhenos~RGenos, family = gaussian)
      fitAlt = glm(ContPhenos~RGenos+Score.gene, family = gaussian, epsilon = 1e-6)	
    }
  } else {
    #if fitting with covariates:
    CovData = allData[[6]]
    CatPhenos = allData[[7]]
    ContPhenos = allData[[8]]
    
    if(Cat == TRUE){
      #fit the null and alternative models
      fitNull = glm(CatPhenos~CovData+RGenos, family = binomial)
      fitAlt = glm(CatPhenos~CovData+RGenos+Score.gene, family = binomial, epsilon = 1e-6)
    } else {
      #continuous outcome fit
      fitNull = glm(ContPhenos~CovData+RGenos, family = gaussian)
      fitAlt = glm(ContPhenos~CovData+RGenos+Score.gene, family = gaussian, epsilon = 1e-6)	
    }
   }
  
  nullValues = summary(fitNull)$coefficients
  altValues = summary(fitAlt)$coefficients
 
  if(fitNull$deviance - fitAlt$deviance >= 0){
    # construct the likelihood ratio statistics
    TL = lrtest(fitNull, fitAlt)
    
    pv = TL[2,5]
    stat = TL[2,4]
    doF = TL[2,1]
  } else {
    pv = 1
    stat = 0
    doF = (dim(RGenos)[2]+1)
  }
  statsAndPValsmat = c(pv, stat)
  output = list(statsAndPValsmat, nullValues, altValues)
  return(output)
}

RunLHT_Cat = function(allData, FitScore, covs = TRUE){
  pv <- rep(0, 3)
  Tall <- matrix(0, 1, 3)
  beta.al <- matrix(0, 1, 2)
  
  #pull elements from allData list
  RGenos = allData[[1]]
  IBS.gene = allData[[2]]
  Incomp.gene = allData[[3]]
  AMS.gene = allData[[4]]
  BinMM.gene = allData[[5]]
  
  if(FitScore == "IBS"){
    Score.gene = IBS.gene
  } else if(FitScore == "Incomp"){
    Score.gene = Incomp.gene
  } else if(FitScore == "AMS"){
    Score.gene = AMS.gene
  } else{
    Score.gene = BinMM.gene
  }
  
  if(covs == FALSE){
    #if no covariates:
    CatPhenos = allData[[6]]
    numCov = 0
	  unpen=c()
  } else {
    #if fitting with covariates:
    CovData = allData[[6]]
    CatPhenos = allData[[7]]
    numCov = dim(CovData)[2] + 1
    #need to include column of 1s for intercept?
    intercept = matrix(1,nrow = numPairs, ncol = 1)
	  unpen=c(1,2,3)
  }
  # define location of zero components
  numSNPs = dim(RGenos)[2]
  N = c(rep(FALSE, numCov+numSNPs), TRUE)
  
  #need to combine the covariates, r genos, and score matrices together
  if(covs == FALSE){
    designMat = cbind(RGenos, Score.gene)
  } else {
    designMat = cbind(intercept, CovData, RGenos, Score.gene)
  }
  #then combine the phenos with the design matrix as a list
  Model = list(X=designMat, Y=CatPhenos)
  
  # estimate the uncontrained estimator
  beta.unre <- cv.SCAD_ADMM_unre(X=Model$X, Y=Model$Y, N=N, beta0=rep(0, dim(Model$X)[2]), err=1e-4, tune="cv", unpen=unpen)
  indice.unre <- beta.unre!=0
  # estimate the constrained estimator (only need this for score test)
  beta.re <- cv.SCAD_ADMM_re(X=Model$X, Y=Model$Y, N=N, beta0=rep(0, dim(Model$X)[2]), err=1e-4, tune="cv", unpen=unpen)
  indice.re <- beta.re!=0
  
  pi.unre <- logit(Model$X%*%beta.unre)
  pi.re <- logit(Model$X%*%beta.re)
  
  # construct the likelihood ratio statistics
  TL <- 2*(sum(log(1+exp(Model$X%*%beta.re))-Model$Y*(Model$X%*%beta.re))-
             sum(log(1+exp(Model$X%*%beta.unre))-Model$Y*(Model$X%*%beta.unre)))
  #if LRT stat is less than 0, there is error, so set to -1
  if(TL < 0){
    TL = -1
  }
  
  # construct the Wald statistics
  #B_0 should give you Omega_a hat
  A = crossprod(Model$X[,indice.unre|N], as.vector(pi.unre*(1-pi.unre))*Model$X[,indice.unre|N])
  if(rcond(A) >= 1e-10){
    B_01 <- solve(A)
    #d_0 should determine which rows to subset
    #gives the first value that is restricted
    d_01 <- sum(which(indice.unre|N)<=length(N))
    #so the B_0 needs to be subset to m rows and columns that are restricted based on H0
    TW <- crossprod(beta.unre[N], solve(B_01[d_01,d_01], beta.unre[N]))
  } else {
    TW = -1
  }
  
  # construct the score statistics
  eps <- Model$Y-pi.re #Y - e(Y)
  #this is the X^T times Y-E(Y)
  Xeps <- crossprod(Model$X[,indice.re|N], eps)
  
  C = crossprod(Model$X[,indice.re|N], as.vector(pi.re*(1-pi.re))*Model$X[,indice.re|N])
  if(rcond(C) >= 1e-10){
    TS <- crossprod(Xeps, solve(C, Xeps))
  } else {
    TS = -1
  }
  
  #determine if null hyp is rejected or not
  #degrees of freedom is equal to the number of restricted values
  doF = sum(N == TRUE)
  
  #set pv to 1 if reject, 0 otherwise
  #LRT
  if (TL>=qchisq(0.95, df = doF)){
    pv[1] <- pv[1]+1/5000
  }
  #Wald
  if (TW>=qchisq(0.95, df = doF)){
    pv[2] <- pv[2]+1/5000
  }
  #Score
  if (TS>=qchisq(0.95, df = doF)){
    pv[3] <- pv[3]+1/5000
  }
  
  Tall[1, ] <- c(TL, TW, TS)
  beta.al[1, ] <- c(sum(beta.re!=0), sum(beta.unre!=0))
  
  statsAndPVals = list(pv=pv, TScores=Tall, beta=beta.al)
  return(statsAndPVals)
}

RunLHT_Cont = function(allData, FitScore, covs = TRUE){
  pv <- rep(0, 3)
  Tall <- matrix(0, 1, 3)
  beta.al <- matrix(0, 1, 2)
  
  #pull elements from allData list
  RGenos = allData[[1]]
  IBS.gene = allData[[2]]
  Incomp.gene = allData[[3]]
  AMS.gene = allData[[4]]
  BinMM.gene = allData[[5]]
  
  if(FitScore == "IBS"){
    Score.gene = IBS.gene
  } else if(FitScore == "Incomp"){
    Score.gene = Incomp.gene
  } else if(FitScore == "AMS"){
    Score.gene = AMS.gene
  } else{
    Score.gene = BinMM.gene
  }
  
  if(covs == FALSE){
    #if no covariates:
    ContPhenos = allData[[7]]
    numCov = 0
	unpen=c()
  } else {
    #if fitting with covariates:
    CovData = allData[[6]]
    ContPhenos = allData[[8]]
    numCov = dim(CovData)[2]
	unpen=c(1,2,3)
  }
  # define location of zero components
  numSNPs = dim(RGenos)[2]
  N = c(rep(FALSE, numCov+1+numSNPs), TRUE)
  #need to include column of 1s for intercept?
  intercept = matrix(1,nrow = numPairs, ncol = 1)
  
  #need to combine the covariates, r genos, and score matrices together
  if(covs == FALSE){
    designMat = cbind(intercept, RGenos, Score.gene)
  } else {
    designMat = cbind(intercept, CovData, RGenos, Score.gene)
  }
  #then combine the phenos with the design matrix as a list
  Model = list(X=designMat, Y=ContPhenos)
  
  # estimate the uncontrained estimator
  beta.unre <- cv.SCAD_ADMM_unre(X=Model$X, Y=Model$Y, N=N, beta0=rep(0, dim(Model$X)[2]), err=1e-4, tune="cv", unpen = unpen)
  indice.unre <- beta.unre!=0
  
  # estimate the constrained estimator
  beta.re <- cv.SCAD_ADMM_re(X=Model$X, Y=Model$Y, N=N, beta0=rep(0, dim(Model$X)[2]), err=1e-4, tune="cv", unpen = unpen)
  indice.re <- beta.re!=0
  
  # estimate the conditional variance
  #    sig2 <- mean((Model$Y-Model$X%*%beta.unre)^2)
  # n = numSims
  # sig2 <- mean((Model$Y-Model$X%*%beta.unre)^2)*n/(n-sum(beta.unre!=0))
  
  # construct the likelihood ratio statistic
  TL <- sum((Model$X%*%beta.re-Model$Y)^2)-sum((Model$X%*%beta.unre-Model$Y)^2)
  if(TL < 0){
    TL = -1
  }
  
  # construct the Wald statistic
  B = crossprod(Model$X[,indice.unre|N], Model$X[,indice.unre|N])
  if(rcond(B) >= 1e-10){
    B_0 <- solve(B)
    d_0 <- sum(which(indice.unre|N)<=length(N))
    TW <- crossprod(beta.unre[N], solve(B_0[d_0,d_0], beta.unre[N]))
  } else {
    TW = -1
  }
  
  # construct the score statistic
  eps <- Model$Y-Model$X%*%beta.re
  Xeps <- crossprod(Model$X[,indice.re|N], eps)
  D = crossprod(Model$X[,indice.re|N], Model$X[,indice.re|N])
  if(rcond(D) >= 1e-10){
    TS <- crossprod(Xeps, solve(D, Xeps))
  } else {
    TS = -1
  }
  
  #determine if null hyp is rejected or not
  #degrees of freedom is equal to the number of restricted values
  doF = sum(N == TRUE)
  
  #set pv to 1 if reject, 0 otherwise
  #LRT
  if (TL>=qchisq(0.95, df = doF)){
    pv[1] <- pv[1]+1/5000
  }
  #Wald
  if (TW>=qchisq(0.95, df = doF)){
    pv[2] <- pv[2]+1/5000
  }
  #Score
  if (TS>=qchisq(0.95, df = doF)){
    pv[3] <- pv[3]+1/5000
  }
  
  Tall[1, ] <- c(TL, TW, TS)
  beta.al[1, ] <- c(sum(beta.re!=0), sum(beta.unre!=0))
  
  statsAndPVals = list(pv=pv, TScores=Tall, beta=beta.al)
  return(statsAndPVals)
}

OrganizeOutput = function(statsAndPVals){
  statsAndPValsAll = list()
  nulValuesAll = list()
  altValuesAll = list()
  
  for(ii in 1:length(statsAndPVals)){
    statsAndPValsAll[[ii]] = statsAndPVals[[ii]][1]
    nulValuesAll[[ii]] =  statsAndPVals[[ii]][2]
    altValuesAll[[ii]] = statsAndPVals[[ii]][3]
  }
  
  listLength = length(statsAndPValsAll)
  statsAndPValsAll.mat = matrix(unlist(statsAndPValsAll),nrow = listLength, ncol = 2, byrow = TRUE)
  dim1 = dim(nulValuesAll[[1]][[1]])
  dim2 = dim(altValuesAll[[1]][[1]])
  rownames1 = c()
  rownames2 = c()
  for(jj in 1:length(statsAndPVals)){
    rownames1 = c(rownames1, rownames(nulValuesAll[[jj]][[1]]))
    rownames2 = c(rownames2, rownames(altValuesAll[[jj]][[1]]))
  }
  nulValuesAll.mat = matrix(unlist(nulValuesAll), ncol = 4, byrow = TRUE, dimnames = list(rownames1, colnames(nulValuesAll[[1]][[1]])))
  altValuesAll.mat = matrix(unlist(altValuesAll), ncol = 4, byrow = TRUE, dimnames = list(rownames2, colnames(altValuesAll[[1]][[1]])))
  
  reorgData = list(statsAndPValsAll.mat,nulValuesAll.mat,altValuesAll.mat)
  return(reorgData)
}

OrganizeOutputPower = function(statsAndPVals){
  n = length(statsAndPVals)
  
  pVals.1 = pVals.2 = pVals.3 = pVals.4 = matrix(nrow = n, ncol = 2)
  NullValues.1 = NullValues.2 = NullValues.3 = NullValues.4 = list()
  altValues.1 = altValues.2 = altValues.3 = altValues.4 = list()
  Score.1 = Score.2 = Score.3 = Score.4 = matrix(nrow = 1, ncol = 1)
  
  for(ll in 1:n){
    pVals.1[ll,] = unlist(statsAndPVals[[ll]][[1]])
    NullValues.1[[ll]] = unlist(statsAndPVals[[ll]][[2]])
    altValues.1[[ll]] = unlist(statsAndPVals[[ll]][[3]])
    Score.1[1,] = unlist(statsAndPVals[[1]][[4]])
    
    pVals.2[ll,] = unlist(statsAndPVals[[ll]][[5]])
    NullValues.2[[ll]] = unlist(statsAndPVals[[ll]][[6]])
    altValues.2[[ll]] = unlist(statsAndPVals[[ll]][[7]])
    Score.2[1,] = unlist(statsAndPVals[[1]][[8]])
    
    pVals.3[ll,] = unlist(statsAndPVals[[ll]][[9]])
    NullValues.3[[ll]] = unlist(statsAndPVals[[ll]][[10]])
    altValues.3[[ll]] = unlist(statsAndPVals[[ll]][[11]])
    Score.3[1,] = unlist(statsAndPVals[[1]][[12]])
    
    pVals.4[ll,] = unlist(statsAndPVals[[ll]][[13]])
    NullValues.4[[ll]] = unlist(statsAndPVals[[ll]][[14]])
    altValues.4[[ll]] = unlist(statsAndPVals[[ll]][[15]])
    Score.4[1,] = unlist(statsAndPVals[[1]][[16]])
  }
  
  pVals = cbind(pVals.1, pVals.2, pVals.3, pVals.4)
  nulValues.1.mat = do.call(rbind,NullValues.1)
  altValues.1.mat = do.call(rbind,altValues.1)
  
  nulValues.2.mat = do.call(rbind,NullValues.2)
  altValues.2.mat = do.call(rbind,altValues.2)
  
  nulValues.3.mat = do.call(rbind,NullValues.3)
  altValues.3.mat = do.call(rbind,altValues.3)
  
  nulValues.4.mat = do.call(rbind,NullValues.4)
  altValues.4.mat = do.call(rbind,altValues.4)
  
  Scores.mat = c(Score.1, Score.2, Score.3, Score.4)
  
  reorgData = list(pVals,nulValues.1.mat,altValues.1.mat,nulValues.2.mat,altValues.2.mat,nulValues.3.mat,altValues.3.mat, nulValues.4.mat,altValues.4.mat,Scores.mat)
  return(reorgData)
}

RunDesparseLasso = function(allData, FitScore, Cat = TRUE, covs = TRUE){
  suppressMessages(require(hdi))
  
  #pull elements from allData list
  RGenos = allData[[1]]
  IBS.gene = allData[[2]]
  Incomp.gene = allData[[3]]
  AMS.gene = allData[[4]]
  BinMM.gene = allData[[5]]
  
  if(FitScore == "IBS"){
    Score.gene = IBS.gene
  } else if(FitScore == "Incomp"){
    Score.gene = Incomp.gene
  } else if(FitScore == "AMS"){
    Score.gene = AMS.gene
  } else{
    Score.gene = BinMM.gene
  }
  
  if(covs == FALSE){
    #if no covariates:
    CatPhenos = allData[[6]]
    ContPhenos = allData[[7]]
    #design mat, no covs
    dMat = cbind(RGenos, Score.gene)
    
    if(Cat == TRUE){
      lassoModel = lasso.proj(x = dMat, y = CatPhenos, family = "binomial", suppress.grouptesting = TRUE)
    } else {
      lassoModel = lasso.proj(x = dMat, y = ContPhenos, suppress.grouptesting = TRUE)
    }
  } else {
    #if fitting with covariates:
    CovData = allData[[6]]
    CatPhenos = allData[[7]]
    ContPhenos = allData[[8]]
    #design mat with covs:
    dMat = cbind(CovData, RGenos, Score.gene)
    
    if(Cat == TRUE){
      lassoModel = lasso.proj(x = dMat, y = CatPhenos, family = "binomial", suppress.grouptesting = TRUE)
    } else {
      lassoModel = lasso.proj(x = dMat, y = ContPhenos, suppress.grouptesting = TRUE)
    }
  }
  scorelocation = length(lassoModel$pval)
  score.pval = lassoModel$pval[scorelocation]
  return(score.pval)
}

RunAISPU = function(allData, FitScore, Cat = TRUE, covs = TRUE, seed){
  #had to change penalty to scad for NAT2 500, see if tlp works in other situations!
  #trying tlp for NAT2 1000 pairs (still doesn't work, sticking with scad)
  
  suppressMessages(require(aispu,lib.loc='/home/vlynn/R/library'))
  
  #pull elements from allData list
  RGenos = allData[[1]]
  IBS.gene = allData[[2]]
  Incomp.gene = allData[[3]]
  AMS.gene = allData[[4]]
  BinMM.gene = allData[[5]]
  
  if(FitScore == "IBS"){
    Score.gene = IBS.gene
  } else if(FitScore == "Incomp"){
    Score.gene = Incomp.gene
  } else if(FitScore == "AMS"){
    Score.gene = AMS.gene
  } else{
    Score.gene = BinMM.gene
  }
  
  numPairs = dim(RGenos)[1]
  
  #center the score variable
  mean.Score.gene = mean(Score.gene)
  Score.gene.centered = Score.gene - mean.Score.gene
  
  #make RGenos cov2 into matrix
  cov2Data = as.matrix(RGenos, nrow = numPairs)
  
  if(covs == FALSE){
    #if no covariates:
    CatPhenos = allData[[6]]
    ContPhenos = allData[[7]]
    
    if(Cat == TRUE){
      model = tryCatch({
        aispu(Y=CatPhenos, X=Score.gene.centered, cov=NULL, cov2=cov2Data, pow = c(1:6, Inf), model = "binomial", penalty = "scad", n.perm = 1000, resample = "boot")
      },
      error = function(error_condition) {
        #change seed and try again
        set.seed(seed+2500)
        aispu(Y=CatPhenos, X=Score.gene.centered, cov=NULL, cov2=cov2Data, pow = c(1:6, Inf), model = "binomial", penalty = "scad", n.perm = 1000, resample = "boot")
      })
    } else {
      model = tryCatch({
        aispu(Y=ContPhenos, X=Score.gene.centered, cov=NULL, cov2=cov2Data, pow = c(1:6, Inf), model = "gaussian", penalty = "scad", n.perm = 1000, resample = "boot")
      },
      error = function(error_condition) {
        #change seed and try again
        set.seed(seed+2500)
        aispu(Y=ContPhenos, X=Score.gene.centered, cov=NULL, cov2=cov2Data, pow = c(1:6, Inf), model = "gaussian", penalty = "scad", n.perm = 1000, resample = "boot")
      })
    }
  } else {
    #if fitting with covariates:
    CovData = allData[[6]]
    CatPhenos = allData[[7]]
    ContPhenos = allData[[8]]
    
    #make sure it is in matrix form
    Covs = as.matrix(CovData, nrow = numPairs)
    
    if(Cat == TRUE){
      model = tryCatch({
        aispu(Y=CatPhenos, X=Score.gene.centered, cov=Covs, cov2=cov2Data, pow = c(1:6, Inf), model = "binomial", n.perm = 1000, penalty = "scad", resample = "boot")
      },
      error = function(error_condition) {
        #change seed and try again
        set.seed(seed+2500)
        aispu(Y=CatPhenos, X=Score.gene.centered, cov=Covs, cov2=cov2Data, pow = c(1:6, Inf), model = "binomial", n.perm = 1000, penalty = "scad", resample = "boot")
      })        
    } else {
      model = tryCatch({
        aispu(Y=ContPhenos, X=Score.gene.centered, cov=Covs, cov2=cov2Data, pow = c(1:6, Inf), model = "gaussian", n.perm = 1000, penalty = "scad", resample = "boot")
      },
      error = function(error_condition) {
        #change seed and try again
        set.seed(seed+2500)
        aispu(Y=ContPhenos, X=Score.gene.centered, cov=Covs, cov2=cov2Data, pow = c(1:6, Inf), model = "gaussian", n.perm = 1000, penalty = "scad", resample = "boot")
      })          
    }
  }
  pValues = c(model$pvs[7], model$pvs[8])
  scores = c(model$Ts[7], model$Ts[8])
  coefs = model$coef
  AllOutput = list(scores, pValues, coefs)
  return(AllOutput)
}

##functions for decorr score test
estVar = function(n, Y, BetaHat, Q){
  XMat = as.matrix(cbind(rep(1,n),Q),nrow=n)
  bHat = as.matrix(BetaHat)
  varValue = (1/n)*sum((Y - XMat%*%bHat)^2)
  return(varValue)
}

calcSHat = function(var, n, Y, gammaHat, X, Z, omegaHat){
  XMat = as.matrix(cbind(rep(1,n),X),nrow=n)
  gHat = as.matrix(gammaHat)
  wHat = as.matrix(omegaHat)
  S = -(1/(var*n))*sum((Y - XMat%*%gHat)*(Z - XMat%*%wHat))
  return(S)
}

calcIHat = function(var,n,Z,omegaHat,X){
  sumZSq = t(Z)%*%Z
  XMat = as.matrix(cbind(rep(1,n),X),nrow=n)
  sumXZ = t(Z)%*%XMat
  wHat = as.matrix(omegaHat)
  Ihat = (1/var)*((1/n)*sumZSq - (1/n)*sumXZ%*%wHat)
  return(Ihat)
}

expit = function(x){
  a = exp(x)
  value = a/(1+a)
  return(value)
}

calcSHatLogit = function(n, Y, gammaHat, X, Z, omegaHat){
  XMat = as.matrix(cbind(rep(1,n),X),nrow=n)
  gHat = as.matrix(gammaHat)
  wHat = as.matrix(omegaHat)
  eTerm = expit(XMat%*%gHat)
  term1 = (Y - eTerm)
  term2 = (Z - XMat%*%wHat)
  S = -(1/n)*sum(term1*term2)
  return(S)
}

calcIHatLogit = function(n,betaHat,Z,omegaHat,X){
  XMat = as.matrix(cbind(rep(1,n),X),nrow=n)
  Q = cbind(XMat, Z)
  bHat = as.matrix(betaHat)
  wHat = as.matrix(omegaHat)
  QBetaHat = Q%*%bHat
  term1 = expit(QBetaHat)*(1/(1+exp(QBetaHat)))
  term3 = (Z - XMat%*%wHat)
  Ihat = (1/n)*sum(term1*Z*term3)
  return(Ihat)
}

RunDecorrScore = function(allData, FitScore, Cat = TRUE, covs = TRUE){
  suppressMessages(require(ncvreg))
  
  #pull elements from allData list
  RGenos = allData[[1]]
  IBS.gene = allData[[2]]
  Incomp.gene = allData[[3]]
  AMS.gene = allData[[4]]
  BinMM.gene = allData[[5]]
  
  if(FitScore == "IBS"){
    Score.gene = IBS.gene
  } else if(FitScore == "Incomp"){
    Score.gene = Incomp.gene
  } else if(FitScore == "AMS"){
    Score.gene = AMS.gene
  } else{
    Score.gene = BinMM.gene
  }
  
  numPairs = dim(RGenos)[1]
  
  if(covs == FALSE){
    #if no covariates:
    CatPhenos = allData[[6]]
    ContPhenos = allData[[7]]
    #design mat, no covs
    dMat = cbind(RGenos, Score.gene)
    
    if(Cat == TRUE){
      betaHats = cv.ncvreg(X=dMat, y = CatPhenos, family = "binomial", penalty = "SCAD")
      betaHatsValues = betaHats$fit$beta[,dim(betaHats$fit$beta)[2]]
      ThetaHat = betaHatsValues[length(betaHatsValues)]
      GammaHat = betaHatsValues[1:(length(betaHatsValues)-1)]
      
      #calc omega hat value
      wHats = cv.ncvreg(X=RGenos, y = Score.gene, family = "gaussian", penalty = "lasso")
      wHatsValues = wHats$fit$beta[,dim(wHats$fit$beta)[2]]
      
      #calculate Shat
      sHat = calcSHatLogit(numPairs, CatPhenos, GammaHat, RGenos, Score.gene, wHatsValues)
      
      #calc Ihat
      iHat = calcIHatLogit(numPairs, betaHatsValues, Score.gene, wHatsValues, RGenos)
      
      #calculate Score Test Stat U
      Uhat = (numPairs)^(1/2)*sHat*(iHat)^(-1/2)
    } else {
      #calc Beta.hat and partition into theta and gamma hat
      betaHats = cv.ncvreg(X=dMat, y = ContPhenos, family = "gaussian", penalty = "SCAD")
      betaHatsValues = betaHats$fit$beta[,dim(betaHats$fit$beta)[2]]
      ThetaHat = betaHatsValues[length(betaHatsValues)]
      GammaHat = betaHatsValues[1:(length(betaHatsValues)-1)]
      
      #calc omega hat value
      wHats = cv.ncvreg(X=RGenos, y = Score.gene, family = "gaussian", penalty = "lasso")
      wHatsValues = wHats$fit$beta[,dim(wHats$fit$beta)[2]]
      
      #estimate variance
      varValue = estVar(numPairs, ContPhenos, betaHatsValues, dMat)
      
      #calc Score
      sHat = calcSHat(varValue, numPairs, ContPhenos, GammaHat, RGenos, Score.gene, wHatsValues)
      
      #calc Ihat for score
      iHat = calcIHat(varValue, numPairs, Score.gene, wHatsValues, RGenos)
      
      #calc score stat
      Uhat = sqrt(numPairs)*sHat*(iHat)^(-1/2)
    }
  } else {
    #if fitting with covariates:
    CovData = allData[[6]]
    CatPhenos = allData[[7]]
    ContPhenos = allData[[8]]
    #design mat with covs:
    dMat = cbind(CovData, RGenos, Score.gene)
    
    if(Cat == TRUE){
      betaHats = cv.ncvreg(X=dMat, y = CatPhenos, family = "binomial", penalty = "SCAD", max.iter=100000)
      betaHatsValues = betaHats$fit$beta[,dim(betaHats$fit$beta)[2]]
      ThetaHat = betaHatsValues[length(betaHatsValues)]
      GammaHat = betaHatsValues[1:(length(betaHatsValues)-1)]
      
      colnames(CovData) = c("Covs1", "Covs2")
      Xmatrix.covs = as.matrix(cbind(CovData, RGenos),nrow = numPairs)
      
      #calc omega hat value
      wHats = cv.ncvreg(X=Xmatrix.covs, y = Score.gene, family = "gaussian", penalty = "lasso", max.iter=100000)
      wHatsValues = wHats$fit$beta[,dim(wHats$fit$beta)[2]]
      
      #calculate Shat
      sHat = calcSHatLogit(numPairs, CatPhenos, GammaHat, Xmatrix.covs, Score.gene, wHatsValues)
      
      #calc Ihat
      iHat = calcIHatLogit(numPairs, betaHatsValues, Score.gene, wHatsValues, Xmatrix.covs)
      
      #calculate Score Test Stat U
      Uhat = (numPairs)^(1/2)*sHat*(iHat)^(-1/2)
    } else {
      #calc Beta.hat and partition into theta and gamma hat
      betaHats = cv.ncvreg(X=dMat, y = ContPhenos, family = "gaussian", penalty = "SCAD", max.iter=100000)
      betaHatsValues = betaHats$fit$beta[,dim(betaHats$fit$beta)[2]]
      ThetaHat = betaHatsValues[length(betaHatsValues)]
      GammaHat = betaHatsValues[1:(length(betaHatsValues)-1)]
      
      colnames(CovData) = c("Covs1", "Covs2")
      Xmatrix.covs = as.matrix(cbind(CovData, RGenos),nrow = numPairs)
      
      #calc omega hat value
      wHats = cv.ncvreg(X=Xmatrix.covs, y = Score.gene, family = "gaussian", penalty = "lasso", max.iter=100000)
      wHatsValues = wHats$fit$beta[,dim(wHats$fit$beta)[2]]
      
      #estimate variance
      varValue = estVar(numPairs, ContPhenos, betaHatsValues, dMat)
      
      #calc Score
      sHat = calcSHat(varValue, numPairs, ContPhenos, GammaHat, Xmatrix.covs, Score.gene, wHatsValues)
      
      #calc Ihat for score
      iHat = calcIHat(varValue, numPairs, Score.gene, wHatsValues, Xmatrix.covs)
      
      #calc score stat
      Uhat = sqrt(numPairs)*sHat*(iHat)^(-1/2)
    }
  }
  p.val = pnorm(abs(Uhat),lower.tail=FALSE)
  score = Uhat
  Output = list(score, p.val)
  return(Output)
}