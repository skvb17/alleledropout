#getting database "genos"
getwd()
setwd("/Users/shantothenel/Desktop/R")
load("/Users/shantothenel/Desktop/R/Rori's Database.Rdata")

#getting allele frequency table
freq = read.table("AFs_AfA.dat.txt") #reads in data
t(freq) #makes a table out of it

#getting random profile for experiment and turning it into an array
ngenos = dim(genos)[1]
sample(ngenos,1) #chose person 30732's profile as our control
controlProfile <- genos[30732,,] #THIS IS OUR PROFILE THAT WE WILL DROP OUT ALLELES FROM
sampleProf <- genos[111425,,] #for testing fxns later


#getting function that counts how many alleles match between two genotype arrays
alleleMatch = function(arr1, arr2) {
  matches = 0
  for(i in 1:13) {
    if((arr1[1,i] == arr2[1,i]) && (arr1[2,i] == arr2[2,i]) || (arr1[1,i] == arr2[2,i] && arr1[2,i] == arr2[1,i]))
      matches = matches + 2
    else if((arr1[1,i] == arr2[1,i]) || (arr1[1,i] == arr2[2,i]) || (arr1[2,i] == arr2[1,i]) || (arr1[2,i] == arr2[2,i]))
      matches = matches + 1
  }
  
  return(matches)
}

alleleMatch(controlProfile,sampleProf)


#getting function that counts the allelic matches between a profile and the first n profiles in a database
directMatch = function(array, database, nProfiles) {
  directmatches = rep(0,27)
  for (i in 1:nProfiles) {
    directmatches0 = alleleMatch(array,database[i,,])
    directmatches[directmatches0 + 1] = directmatches[directmatches0 + 1] + 1
  }
  
  return(directmatches)
}

#entering the control genotype into function directMatch
directMatch(controlProfile, genos, ngenos)

#yielded a full match at 26 out of 26 alleles matching; no to low chance of false positives



########## SIMULATING DROPOUT AT 1 ALLELE ##########

##### first trial #####

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#function for ONE allele dropout at a random loci SOFTCODE
  #input: selected genotype (profile)
  #purpose: changes matrix so that it changes chosen allele to 0 effectively dropping an allele out
  #output: genotype with one random allele dropped out

dropout1 = function(profile) {
  locidrop <- sample(13,1,F)
  alledrop <- sample(2,1,T) 
  profile[alledrop, locidrop] <- 0
  
  return(profile)
}

controlProf01a <- dropout1(controlProfile) #works!!

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#entering the FIRST dropout genotype with ONE allele into function directMatch
directMatch(controlProf01a, genos, ngenos)

#yielded ONE partial match at 25 out of 26 alleles matching; low chance of false positive
#yielded TWO partial matches at 18 out of 26 alleles matching
#yielded FIFTEEN partial matches at 17 out of 26 alleles matching



##### second trial #####
controlProf01b <- dropout1(controlProfile) #works!!

#entering the SECOND dropout genotype with ONE allele into function directMatch
directMatch(controlProf01b, genos, ngenos)

#yielded ONE partial match at 25 out of 26 alleles matching; low chance of false positive
#yielded FOUR partial matches at 18 out of 26 alleles matching and THIRTEEN at 17 out of 26 alleles matching



##### third trial #####
controlProf01c <- dropout1(controlProfile) #works!!

#entering the THIRD dropout genotype with ONE allele into function directMatch
directMatch(controlProf01c, genos, ngenos)

#yielded ONE partial match at 25 out of 26 alleles matching; low chance of false positive
#yielded ONE partial match at 18 out of 26 alleles matching and FIFTEEN at 17 out of 26 alleles matching



##### fourth trial #####
controlProf01d <- dropout1(controlProfile) #works!!

#entering the FOURTH dropout genotype with ONE allele into function directMatch
directMatch(controlProf01d, genos, ngenos)

#yielded ONE partial match at 25 out of 26 alleles matching; low chance of false positive
#yielded ONE partial match at 18 out of 16 alleles matching and ELEVEN at 17 out of 26 alleles matching



##### fifth trial #####
controlProf01e <- dropout1(controlProfile) #works!!

#entering the FIFTH dropout genotype with ONE allele into function directMatch
directMatch(controlProf01e, genos, ngenos)

#yielded ONE partial match at 25 out of 26 alleles matching; low chance of false positive
#yielded ONE partial match at 18 out of 26 alleles matching and TEN at 17 our of 26 alleles matching



########## SIMULATING DROPOUT AT 5 ALLELES ##########

##### first trial #####

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#function for FIVE allele dropouts at a random loci SOFTCODE
#input: selected genotype (profile)
#purpose: changes matrix so that it changes chosen allele to 0 effectively dropping an allele out
#output: genotype with five random alleles dropped out

dropout5 = function(profile) {
  locidrop <- sample(13,5,F)
  alledrop <- sample(2,5,T) 
  for (i in 1:5)
    if(profile[alledrop[i], locidrop[i]] <- 0) {
      print(profile)
    }

  return(profile)
}

controlProf05a <- dropout5(controlProfile) #works!!

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#entering the FIRST dropout genotype with FIVE alleles into function directMatch
directMatch(controlProf05a, genos, ngenos)

#yielded ONE partial match at 21 out of 26 alleles matching; low chance of false positive
#yielded TWO partial matches at 16 out of 26 alleles matching



##### second trial #####
controlProf05b <- dropout5(controlProfile) #works!!

#entering the SECOND dropout genotype with FIVE alleles into function directMatch
directMatch(controlProf05b, genos, ngenos)

#yielded ONE partial match at 21 out of 26 alleles matching; low chance of false positive
#yielded TEN partial matches at 15 out of 26 alleles matching



##### third trial #####
controlProf05c <- dropout5(controlProfile)

#entering the THIRD dropout genotype with FIVE alleles into function directMatch
directMatch(controlProf05c, genos, ngenos)

#yielded ONE partial match at 21 out of 26 alleles matching; low chance of false positive
#yielded SIX partial match at 16 out of 26 alleles matching



##### fourth trial #####
controlProf05d <- dropout5(controlProfile)

#entering the FOURTH dropout genotype with FIVE alleles into function directMatch
directMatch(controlProf05d, genos, ngenos)

#yielded ONE partial match at 21 out of 26 alleles matching; low chance of false positive
#yielded FIVE partial match at 15 out of 26 alleles matching



##### fifth trial #####
controlProf05e <- dropout5(controlProfile)

#entering the FIFTH dropout genotype with FIVE alleles into function directMatch
directMatch(controlProf05e, genos, ngenos)

#yielded ONE partial match at 21 out of 26 alleles matching; low chance of false positive
#yielded FIFTEEN partial match at 15 out of 26 alleles matching



########## SIMULATING DROPOUT AT 7 ALLELES ##########

##### first trial #####

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#function for SEVEN allele dropouts at a random loci SOFTCODE
#input: selected genotype (profile)
#purpose: changes matrix so that it changes chosen allele to 0 effectively dropping an allele out
#output: genotype with seven random alleles dropped out

dropout7 = function(profile) {
  locidrop <- sample(13,7,F)
  alledrop <- sample(2,7,T) 
  for (i in 1:7)
    if(profile[alledrop[i], locidrop[i]] <- 0) {
      print(profile)
    }
  
  return(profile)
}

controlProf07a <- dropout7(controlProfile) #works!!

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#entering the FIRST dropout genotype with SEVEN alleles into function directMatch
directMatch(controlProf07a, genos, ngenos)

#yielded ONE partial match at 19 out of 26 alleles matching; low chance of false positive
#yielded SIX partial matches at 14 out of 26 alleles matching



##### second trial #####
controlProf07b <- dropout7(controlProfile)

#entering the SECOND dropout genotype with SEVEN alleles into function directMatch
directMatch(controlProf07b, genos, ngenos)

#yielded ONE partial match at 19 out of 26 alleles matching; low chance of false positive
#yielded ONE partial match at 15 out of 26 alleles matching



##### third trial #####
controlProf07c <- dropout7(controlProfile)

#entering the THIRD dropout genotype with SEVEN alleles into function directMatch
directMatch(controlProf07c, genos, ngenos)

#yielded ONE partial match at 19 out of 26 alleles matching; low chance of false positive
#yielded NINE partial matches at 14 out of 26 alleles matching



##### fourth trial #####
controlProf07d <- dropout7(controlProfile)

#entering the FOURTH dropout genotype with SEVEN alleles into function directMatch
directMatch(controlProf07d, genos, ngenos)

#yielded ONE partial match at 19 out of 26 alleles matching; low chance of false positive
#yielded ONE partial matches at 16 out of 26 alleles matching



##### fifth trial #####
controlProf07e <- dropout7(controlProfile)

#entering the FIFTH dropout genotype with SEVEN alleles into function directMatch
directMatch(controlProf07e, genos, ngenos)

#yielded ONE partial match at 19 out of 26 alleles matching; low chance of false positive
#yielded FIFTEEN partial matches at 14 out of 26 alleles matching



########## SIMULATING DROPOUT AT 10 ALLELES ##########

##### first trial #####

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#function for TEN allele dropouts at a random loci SOFTCODE
#input: selected genotype (profile)
#purpose: changes matrix so that it changes chosen allele to 0 effectively dropping an allele out
#output: genotype with ten random alleles dropped out

dropout10 = function(profile) {
  locidrop <- sample(13,10,F)
  alledrop <- sample(2,10,T) 
  for (i in 1:10)
    if(profile[alledrop[i], locidrop[i]] <- 0) {
      print(profile)
    }
  
  return(profile)
}

controlProf10a <- dropout10(controlProfile) #works!!

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#entering the FIRST dropout genotype with TEN alleles into function directMatch
directMatch(controlProf10a, genos, ngenos)

#yielded ONE partial match at 16 out of 26 alleles matching; low chance of false positive
#yielded TWO partial matches at 14 out of 26 alleles matching



##### second trial #####
controlProf10b <- dropout10(controlProfile)

#entering the SECOND dropout genotype with TEN alleles into function directMatch
directMatch(controlProf10b, genos, ngenos)

#yielded ONE partial match at 16 out of 26 alleles matching; low chance of false positive
#yielded ONE partial matches at 13 out of 26 alleles matching



##### third trial #####
controlProf10c <- dropout10(controlProfile)

#entering the THIRD dropout genotype with TEN alleles into function directMatch
directMatch(controlProf10c, genos, ngenos)

#yielded ONE partial match at 16 out of 26 alleles matching; low chance of false positive
#yielded ONE partial matches at 13 out of 26 alleles matching



##### fourth trial #####
controlProf10d <- dropout10(controlProfile)

#entering the FOURTH dropout genotype with TEN alleles into function directMatch
directMatch(controlProf10d, genos, ngenos)

#yielded ONE partial match at 16 out of 26 alleles matching; low chance of false positive
#yielded TWO partial matches at 14 out of 26 alleles matching



##### fifth trial #####
controlProf10e <- dropout10(controlProfile)

#entering the FIFTH dropout genotype with TEN alleles into function directMatch
directMatch(controlProf10e, genos, ngenos)

#yielded ONE partial match at 16 out of 26 alleles matching; low chance of false positive
#yielded ONE partial matches at 14 out of 26 alleles matching



########## SIMULATING DROPOUT AT 13 ALLELES ##########

##### first trial #####

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


#function for THIRTEEN allele dropouts at a random loci SOFTCODE
#input: selected genotype (profile)
#purpose: changes matrix so that it changes chosen allele to 0 effectively dropping an allele out
#output: genotype with thirteen random alleles dropped out

dropout13 = function(profile) {
  locidrop <- sample(13,13,F)
  alledrop <- sample(2,13,T) 
  for (i in 1:13)
    if(profile[alledrop[i], locidrop[i]] <- 0) {
      print(profile)
    }
  
  return(profile)
}

controlProf13a <- dropout13(controlProfile)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#entering the FIRST dropout genotype with THIRTEEN alleles into function directMatch
directMatch(controlProf13a, genos, ngenos)

#yielded ONE partial match at 13 out of 26 alleles matching; low to medium chance of false positive
#yielded THIRTY-SEVEN partial matches at 11 out of 26 alleles matching



##### second trial #####
controlProf13b <- dropout13(controlProfile)

#entering the SECOND dropout genotype with THIRTEEN alleles into function directMatch
directMatch(controlProf13b, genos, ngenos)

#yielded ONE partial match at 13 out of 26 alleles matching; low to medium chance of false positive
#yielded FIFTEEN partial matches at 11 out of 26 alleles matching



##### third trial #####
controlProf13c <- dropout13(controlProfile)

#entering the THIRD dropout genotype with THIRTEEN alleles into function directMatch
directMatch(controlProf13c, genos, ngenos)

#yielded ONE partial match at both 13 out of 26 alleles matching; low to medium chance of false positive
#yielded TEN partial matches at 11 out of 26 alleles matching



##### fourth trial #####
controlProf13d <- dropout13(controlProfile)

#entering the FOURTH dropout genotype with THIRTEEN alleles into function directMatch
directMatch(controlProf13d, genos, ngenos)

#yielded ONE partial match at both 13 out of 26 alleles matching; low to medium chance of false positive
#yielded EIGHT partial matches at 11 out of 26 alleles matching



##### fifth trial #####
controlProf13e <- dropout13(controlProfile)

#entering the FIFTH dropout genotype with THIRTEEN alleles into function directMatch
directMatch(controlProf13e, genos, ngenos)

#yielded ONE partial match at 13 out of 26 alleles matching; low to medium chance of false positive
#yielded ONE partial matches at 12 out of 26 alleles matching



########## CALLING FORENSIM ##########
install.packages("forensim")
install.packages("tcltk2")
install.packages("tkrplot")
library(forensim)
?forensim

##### LR for dropoutProfile1a #####
LRdropoutProfile1a = rep(0,13)
for(i in 1:13) {
  LRdropoutProfile1a[i] = LR(Repliste=dropoutProfile1a[,i], Tp=controlProfile[,i], Td=0, Vp=0, Vd=controlProfile[,i], xp=0, xd=1, theta=0.01, prDHet=0.0384615385, prDHom=0, prC=0, freq=strusa$tab$Afri[[i]])$LR
}
LRdropoutProfile1a
data(strusa)

prod(LRdropoutProfile1a) #2.382209e+16
strusa



##### LR for dropoutProfile5a #####
LRdropoutProfile5a = rep(0,13)
for(i in 1:13) {
  LRdropoutProfile5a[i] = LR(Repliste=dropoutProfile5a[,i], Tp=controlProfile[,i], Td=0, Vp=0, Vd=controlProfile[,i], xp=0, xd=1, theta=0.01, prDHet=0.1923076923, prDHom=0, prC=0, freq=strusa$tab$Afri[[i]])$LR
}
LRdropoutProfile5a

prod(LRdropoutProfile5a) #461499282309



##### LR for dropoutProfile7a #####
LRdropoutProfile7a = rep(0,13)
for(i in 1:13) {
  LRdropoutProfile7a[i] = LR(Repliste=dropoutProfile7a[,i], Tp=controlProfile[,i], Td=0, Vp=0, Vd=controlProfile[,i], xp=0, xd=1, theta=0.01, prDHet=0.3846153846, prDHom=0, prC=0, freq=strusa$tab$Afri[[i]])$LR
}
LRdropoutProfile7a

prod(LRdropoutProfile7a) #6703620614



##### LR for dropoutProfile10a #####
LRdropoutProfile10a = rep(0,13)
for(i in 1:13) {
  LRdropoutProfile10a[i] = LR(Repliste=dropoutProfile10a[,i], Tp=controlProfile[,i], Td=0, Vp=0, Vd=controlProfile[,i], xp=0, xd=1, theta=0.01, prDHet=0.38461538, prDHom=0, prC=0, freq=strusa$tab$Afri[[i]])$LR
}
LRdropoutProfile10a

prod(LRdropoutProfile10a) #16681695



##### LR for dropoutProfile13a #####
LRdropoutProfile13a = rep(0,13)
for(i in 1:13) {
  LRdropoutProfile13a[i] = LR(Repliste=dropoutProfile13a[,i], Tp=controlProfile[,i], Td=0, Vp=0, Vd=controlProfile[,i], xp=0, xd=1, theta=0.01, prDHet=0.5, prDHom=0, prC=0, freq=strusa$tab$Afri[[i]])$LR
}
LRdropoutProfile13a

prod(LRdropoutProfile13a) #154769.8

LRdropoutVec = c(2.382209e+16,461499282309,6703620614,16681695,154769.8)
plot(log(LRdropoutVec))


















