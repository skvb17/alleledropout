#getting database "genos"
getwd()
setwd("~/Desktop/Shannel-project")
load("~/Desktop/Shannel-project/Rori's Database.Rdata")

#getting allele frequency table
freq = read.table("AFs_AfA.dat.txt") #reads in data
t(freq)

#getting random profile for experiment and turning it into an array
dim(genos)  # 200K profiles 2 allelles  13 loci
dim(genos)[1] # 200k
ngenos = dim(genos)[1]
sample(ngenos,1) -> PROFILE  # choosing # of RANDOM profile(s)
PROFILE
genos[PROFILE,,] -> CONTROL # this require a manual insert of the choosen profile # 
CONTROL

# NOT SURE IF WE NEED THIS ANYMORE
#controlVec = c(9,14,23,25,6,8,8,10,17,15,15,17,12,11,12,11,15,14,13,12,11,12,14,18,39.0,33.2)
#controlProfile = array(controlVec, dim=c(2,13))
#controlProfile



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
directMatch(CONTROL, genos, ngenos) 

# [1]     0     1    12   143   717  2718  8049 17631 29812 38103 38716 30782 19032  9270  3592  1138   238    41     4     0     0     0     0     0     0     0     1
# THIS RETURN 0-27 NUMBERS, IT REP HOW MANY MATCHES TO EACH NUMBER OF 0-26 ALLE.


#yielded a full match at 26 out of 26 alleles matching; no to low chance of false positives



########## SIMULATING DROPOUT AT 1 ALLELE ##########

##### first trial #####
sample(13,1) -> LOCIDROP
LOCIDROP 
#12th loci
sample(2,1) ->ALLEDROP
ALLEDROP
#1st allele

ngenos = dim(genos)[1]
sample(ngenos,1) -> PROFILE  # choosing # of RANDOM profile(s)
PROFILE
genos[PROFILE,,] -> CONTROL # this require a manual insert of the choosen profile # 
CONTROL

#CONTROL[ALLEDROP,LOCIDROP]
#gsub(CONTROL[ALLEDROP,LOCIDROP], 0, CONTROL) -> CONTROL

CONTROL[ALLEDROP,LOCIDROP] <-0
CONTROL

#entering the control genotype into function directMatch
directMatch(CONTROL, genos, ngenos) 




###  Completed - 12/16/2017 
## 13 alleles drop with for loops

ngenos = dim(genos)[1]
sample(ngenos,1) -> PROFILE  # choosing # of RANDOM profile(s)
PROFILE
genos[PROFILE,,] -> CONTROL # this require a manual insert of the choosen profile # 
CONTROL

DROP13 <- matrix(0, ncol = 27, nrow = 13)
DROP13 <- data.frame(DROP13)
rownames(DROP13) <- c("Drop1","Drop2", "Drop3", "Drop4","Drop5", "Drop6", "Drop7","Drop8", "Drop9", "Drop10","Drop11","Drop12", "Drop13")


result =NULL
directmatchResult=NULL

for (j in c(1:13)){ 
    dropout = function(DDprofile) {
                                    LOCIDROP <- sample(13,j,F)
                                    ALLEDROP <- sample(2,j,T) 
                                    for (i in c(1:j))
                                    if(DDprofile[ALLEDROP[i], LOCIDROP[i]] <- 0) {
                                                                                  print(DDprofile)
                                                                                  }
                                    return(DDprofile)
                                    }
                  result = dropout(CONTROL) 
                  print(result)
                  directmatchResult = unlist(directMatch(result, genos, ngenos)) -> DROP13[j,]
                  #print(unlist(directMatch(result, genos, ngenos)))
                  }


#### trying to run it 20 times for DROP 1 allele ONLY


DROP1 <- matrix(0, ncol = 27, nrow = 10)
DROP1 <- data.frame(DROP1)

result =NULL
directmatchResult=NULL


dropout7 = function(profile) {
  locidrop <- sample(13,7,F)
  alledrop <- sample(2,7,T) 
  for (i in 1:7)
    if(profile[alledrop[i], locidrop[i]] <- 0) {
      print(profile)
    }
  
  return(profile)
}











#genotype controlProfile except with the 12th locus, 1st allele dropped out
dropoutVec1a = c(9,14,23,25,6,8,8,10,17,15,15,17,12,11,12,11,15,14,13,12,11,12,0,18,39.0,33.2)
dropoutProfile1a = array(dropoutVec1a, dim=c(2,13))
dropoutProfile1a

#entering the FIRST dropout genotype with ONE allele into function directMatch
directMatch(dropoutProfile1a, genos, ngenos)

#yielded ONE partial match at 25 out of 26 alleles matching; low chance of false positive
#yielded ONE partial match at both 19 out of 26 and 18 out of 26 alleles matching
#yielded EIGHTEEN partial matches at 17 out of 26 alleles matching

oneA = c(0,0,22,181,1032,4122,11304,23025,35766,41329,36895,25248,13242,5602,1679,433,99,18,1,1,0,0,0,0,0,1,0)
plot(oneA, main = "Trial One's Dropout at One Allele", xlab="Profile Matches", ylab="Alleles", col="blue")
abline(v = 25, col = "red") #indicates number of alleles where the partial match was highest

pdf(file="Trial One's Dropout at One Allele.pdf")
plot(oneA, main = "Trial One's Dropout at One Allele", xlab="Profile Matches", ylab="Alleles", col="blue")
abline(v = 25, col = "red") #indicates number of alleles where the partial match was highest
dev.off()



##### second trial #####
sample(13,1)
#5th loci
sample(2,1)
#1st allele

#genotype controlProfile except with the 5th locus, 1st allele dropped out
dropoutVec1b = c(9,14,23,25,6,8,8,10,0,15,15,17,12,11,12,11,15,14,13,12,11,12,14,18,39.0,33.2)
dropoutProfile1b = array(dropoutVec1b, dim=c(2,13))
dropoutProfile1b

#entering the SECOND dropout genotype with ONE allele into function directMatch
directMatch(dropoutProfile1b, genos, ngenos)

#yielded ONE partial match at 25 out of 26 alleles matching; low chance of false positive
#yielded TWO partial matches at 18 out of 26 alleles matching and ELEVEN at 17 out of 26 alleles matching

oneB = c(0,3,29,305,1632,5722,14516,27523,38893,41083,33421,21243,10296,3842,1135,286,57,11,2,0,0,0,0,0,0,1,0)
plot(oneB, main = "Trial Two's Dropout at One Allele", xlab="Profile Matches", ylab="Alleles", col="blue")
abline(v = 25, col = "red") #indicates number of alleles where the partial match was highest

pdf(file="Trial Two's Dropout at One Allele.pdf")
plot(oneB, main = "Trial Two's Dropout at One Allele", xlab="Profile Matches", ylab="Alleles", col="blue")
abline(v = 25, col = "red") #indicates number of alleles where the partial match was highest
dev.off()



##### third trial #####
sample(13,1)
#8th loci
sample(2,1)
#2nd allele

#genotype controlProfile except with the 8th locus, 2nd allele dropped out
dropoutVec1c = c(9,14,23,25,6,8,8,10,17,15,15,17,12,11,12,0,15,14,13,12,11,12,14,18,39.0,33.2)
dropoutProfile1c = array(dropoutVec1c, dim=c(2,13))
dropoutProfile1c

#entering the THIRD dropout genotype with ONE allele into function directMatch
directMatch(dropoutProfile1c, genos, ngenos)

#yielded ONE partial match at 25 out of 26 alleles matching; low chance of false positive
#yielded TWO partial matches at 18 out of 26 alleles matching and NINE at 17 out of 26 alleles matching

oneC = c(0,0,31,275,1407,5349,13994,26597,38864,41396,34303,21629,10642,3981,1183,268,69,9,2,0,0,0,0,0,0,1,0)
plot(oneC, main = "Trial Three's Dropout at One Allele", xlab="Alleles", ylab="Profile Matches", col="blue")
abline(v = 25, col = "red") #indicates number of alleles where the partial match was highest

pdf(file="Trial Three's Dropout at One Allele.pdf")
plot(oneC, main = "Trial Three's Dropout at One Allele", xlab="Alleles", ylab="Profile Matches", col="blue")
abline(v = 25, col = "red") #indicates number of alleles where the partial match was highest
dev.off()



##### fourth trial #####
sample(13,1)
#4th loci
sample(2,1)
#1st allele

#genotype controlProfile except with the 4th locus, 1st allele dropped out
dropoutVec1d = c(9,14,23,25,6,8,0,10,17,15,15,17,12,11,12,11,15,14,13,12,11,12,14,18,39.0,33.2)
dropoutProfile1d = array(dropoutVec1d, dim=c(2,13))
dropoutProfile1d

#entering the FOURTH dropout genotype with ONE allele into function directMatch
directMatch(dropoutProfile1d, genos, ngenos)

#yielded ONE partial match at 25 out of 26 alleles matching; low chance of false positive
#yielded ONE partial match at 18 out of 16 alleles matching and SEVEN at 17 out of 26 alleles matching

oneD = c(0,2,54,421,1952,6907,16723,30199,40018,40516,31554,18609,8747,3139,885,225,40,7,1,0,0,0,0,0,0,1,0)
plot(oneD, main = "Trial Four's Dropout at One Allele", xlab="Alleles", ylab="Profile Matches", col="blue")
abline(v = 25, col = "red") #indicates number of alleles where the partial match was highest

pdf(file="Trial Four's Dropout at One Allele.pdf")
plot(oneD, main = "Trial Four's Dropout at One Allele", xlab="Alleles", ylab="Profile Matches", col="blue")
abline(v = 25, col = "red") #indicates number of alleles where the partial match was highest
dev.off()



##### fifth trial #####
sample(13,1)
#7th loci
sample(2,1)
#2nd allele

#genotype controlProfile except with the 7th locus, 2nd allele dropped out
dropoutVec1e = c(9,14,23,25,6,8,8,10,17,15,15,17,12,0,12,11,15,14,13,12,11,12,14,18,39.0,33.2)
dropoutProfile1e = array(dropoutVec1e, dim=c(2,13))
dropoutProfile1e

#entering the FIFTH dropout genotype with ONE allele into function directMatch
directMatch(dropoutProfile1e, genos, ngenos)

#yielded ONE partial match at 25 out of 26 alleles matching; low chance of false positive
#yielded ONE partial match at 18 out of 26 alleles matching

oneE = c(0,1,44,343,1641,5897,14781,27369,38561,41099,33821,21003,10179,3855,1080,253,57,14,1,0,0,0,0,0,0,1,0)
plot(oneE, main = "Trial Five's Dropout at One Allele", xlab="Alleles", ylab="Profile Matches", col="blue")
abline(v = 25, col = "red") #indicates number of alleles where the partial match was highest

pdf(file="Trial Five's Dropout at One Allele.pdf")
plot(oneE, main = "Trial Five's Dropout at One Allele", xlab="Alleles", ylab="Profile Matches", col="blue")
abline(v = 25, col = "red") #indicates number of alleles where the partial match was highest
dev.off()



########## SIMULATING DROPOUT AT 5 ALLELES ##########

##### first trial #####
sample(13,5,F)
#13th, 10th, 11th, 7th, and 1st loci respectively
sample(2,5,T)
#2nd, 1st, 2nd, 1st, and 1st alleles respectively

#genotype controlProfile except with the corresponding loci and alleles dropped out
dropoutVec5a = c(0,14,23,25,6,8,8,10,17,15,15,17,0,11,12,11,15,14,0,12,11,0,14,18,39,0)
dropoutProfile5a = array(dropoutVec5a, dim=c(2,13))
dropoutProfile5a

#entering the FIRST dropout genotype with FIVE alleles into function directMatch
directMatch(dropoutProfile5a, genos, ngenos)

#yielded ONE partial match at 21 out of 26 alleles matching; low chance of false positive
#yielded THREE partial matches at 16 out of 26 alleles matching

fiveA = c(1,27,237,1408,5252,14492,28244,40547,42578,34374,20094,8905,3007,695,118,17,3,0,0,0,0,1,0,0,0,0,0)
plot(fiveA, main = "Trial One's Dropout at Five Alleles", xlab="Alleles", ylab="Profile Matches", col="purple")
abline(v = 21, col = "gold") #indicates number of alleles where the partial match was highest

pdf(file="Trial One's Dropout at Five Alleles.pdf")
plot(fiveA, main = "Trial One's Dropout at Five Alleles", xlab="Alleles", ylab="Profile Matches", col="purple")
abline(v = 21, col = "gold") #indicates number of alleles where the partial match was highest
dev.off()



##### second trial #####
sample(13,5,F)
#6th, 9th, 8th, 5th, and 2nd loci respectively
sample(2,5,T)
#2nd, 1st, 2nd, 1st, and 2nd alleles respectively

#genotype controlProfile except with the corresponding loci and alleles dropped out
dropoutVec5b = c(9,14,23,0,6,8,8,10,0,15,0,15,12,11,0,0,14,14,13,12,11,12,14,18,39.0,33.2)
dropoutProfile5b = array(dropoutVec5b, dim=c(2,13))
dropoutProfile5b

#entering the SECOND dropout genotype with FIVE alleles into function directMatch
directMatch(dropoutProfile5b, genos, ngenos)

#yielded ONE partial match at 21 out of 26 alleles matching; low chance of false positive
#yielded SIX partial matches at 15 out of 26 alleles matching

fiveB = c(2,29,357,2083,8005,19578,34851,44153,40658,28174,14649,5466,1597,335,56,6,0,0,0,0,0,1,0,0,0,0,0)
plot(fiveB, main = "Trial Two's Dropout at Five Alleles", xlab="Alleles", ylab="Profile Matches", col="purple")
abline(v = 21, col = "gold") #indicates number of alleles where the partial match was highest

pdf(file="Trial Two's Dropout at Five Alleles.pdf")
plot(fiveB, main = "Trial Two's Dropout at Five Alleles", xlab="Alleles", ylab="Profile Matches", col="purple")
abline(v = 21, col = "gold") #indicates number of alleles where the partial match was highest
dev.off()



##### third trial #####
sample(13,5,F)
#10th, 4th, 9th, 12th, and 11th loci respectively
sample(2,5,T)
#1st, 2nd, 1st, 1st, and 1st alleles respectively

#genotype controlProfile except with the corresponding loci and alleles dropped out
dropoutVec5c = c(9,14,23,25,6,8,8,0,17,15,15,17,12,11,12,11,0,14,0,12,0,12,0,18,39.0,33.2)
dropoutProfile5c = array(dropoutVec5c, dim=c(2,13))
dropoutProfile5c

#entering the THIRD dropout genotype with FIVE alleles into function directMatch
directMatch(dropoutProfile5c, genos, ngenos)

#yielded ONE partial match at 21 out of 26 alleles matching; low chance of false positive
#yielded ONE partial match at 16 out of 26 alleles matching

fiveC = c(0,22,243,1562,5848,15772,30044,41570,42533,32609,18723,7933,2482,532,104,21,1,0,0,0,0,1,0,0,0,0,0)
plot(fiveC, main = "Trial Three's Dropout at Five Alleles", xlab="Alleles", ylab="Profile Matches", col="purple")
abline(v = 21, col = "gold") #indicates number of alleles where the partial match was highest

pdf(file="Trial Three's Dropout at Five Alleles.pdf")
plot(fiveC, main = "Trial Three's Dropout at Five Alleles", xlab="Alleles", ylab="Profile Matches", col="purple")
abline(v = 21, col = "gold") #indicates number of alleles where the partial match was highest
dev.off()



##### fourth trial #####
sample(13,5,F)
#10th, 4th, 9th, 12th, and 11th loci respectively
sample(2,5,T)
#1st, 2nd, 1st, 1st, and 1st alleles respectively

#genotype controlProfile except with the corresponding loci and alleles dropped out
dropoutVec5d = c(9,14,23,25,6,8,8,0,17,15,15,17,12,11,12,11,0,14,0,12,0,12,0,18,39.0,33.2)
dropoutProfile5d = array(dropoutVec5d, dim=c(2,13))
dropoutProfile5d

#entering the FOURTH dropout genotype with FIVE alleles into function directMatch
directMatch(dropoutProfile5d, genos, ngenos)

#yielded ONE partial match at 21 out of 26 alleles matching; low chance of false positive
#yielded ONE partial match at 16 out of 26 alleles matching

fiveD = c(0,22,243,1562,5848,15772,30044,41570,42533,32609,18723,7933,2482,532,104,21,1,0,0,0,0,1,0,0,0,0,0)
plot(fiveD, main = "Trial Four's Dropout at Five Alleles", xlab="Alleles", ylab="Profile Matches", col="purple")
abline(v = 21, col = "gold") #indicates number of alleles where the partial match was highest

pdf(file="Trial Four's Dropout at Five Alleles.pdf")
plot(fiveD, main = "Trial Four's Dropout at Five Alleles", xlab="Alleles", ylab="Profile Matches", col="purple")
abline(v = 21, col = "gold") #indicates number of alleles where the partial match was highest
dev.off()



##### fifth trial #####
sample(13,5,F)
#2nd, 12th, 6th, 8th, and 7th loci respectively
sample(2,5,T)
#1st, 1st, 1st, 2nd, and 1st alleles respectively

#genotype controlProfile except with the corresponding loci and alleles dropped out
dropoutVec5e = c(9,14,0,25,6,8,8,10,17,15,0,17,0,11,12,0,15,14,13,12,11,12,0,18,39.0,33.2)
dropoutProfile5e = array(dropoutVec5e, dim=c(2,13))
dropoutProfile5e

#entering the FIFTH dropout genotype with FIVE alleles into function directMatch
directMatch(dropoutProfile5e, genos, ngenos)

#yielded ONE partial match at 21 out of 26 alleles matching; low chance of false positive
#yielded ONE partial match at both 17 out of 26 and 16 out of 26 alleles matching

fiveE = c(0,47,548,2763,9643,22442,37075,44415,39253,25233,12387,4545,1314,283,44,5,1,1,0,0,0,1,0,0,0,0,0)
plot(fiveE, main = "Trial Five's Dropout at Five Alleles", xlab="Alleles", ylab="Profile Matches", col="purple")
abline(v = 21, col = "gold") #indicates number of alleles where the partial match was highest

pdf(file="Trial Five's Dropout at Five Alleles.pdf")
plot(fiveE, main = "Trial Five's Dropout at Five Alleles", xlab="Alleles", ylab="Profile Matches", col="purple")
abline(v = 21, col = "gold") #indicates number of alleles where the partial match was highest
dev.off()



########## SIMULATING DROPOUT AT 7 ALLELES ##########

##### first trial #####
sample(13,7,F)
#5th, 1st, 12th, 2nd, 11th, 7th, and 4th loci respectively
sample(2,7,T)
#2nd, 1st, 2nd, 1st, 2nd, 1st, and 1st alleles respectively

#genotype controlProfile except with the corresponding loci and alleles dropped out
dropoutVec7a = c(0,14,0,25,6,8,0,10,17,0,15,17,0,11,12,11,15,14,13,12,11,0,14,0,39.0,33.2)
dropoutProfile7a = array(dropoutVec7a, dim=c(2,13))
dropoutProfile7a

#entering the FIRST dropout genotype with SEVEN alleles into function directMatch
directMatch(dropoutProfile7a, genos, ngenos)

#yielded ONE partial match at 19 out of 26 alleles matching; low chance of false positive
#yielded TWO partial matches at 15 out of 26 alleles matching
#yielded THREE partial matches at 14 out of 26 alleles

sevenA = c(8,119,1126,5212,15707,32113,44685,44488,31689,16466,6197,1812,321,51,3,2,0,0,0,1,0,0,0,0,0,0,0)
plot(sevenA, main = "Trial One's Dropout at Seven Alleles", xlab="Alleles", ylab="Profile Matches", col="green")
abline(v = 19, col = "sky blue") #indicates number of alleles where the partial match was highest

pdf(file="Trial One's Dropout at Seven Alleles.pdf")
plot(sevenA, main = "Trial One's Dropout at Seven Alleles", xlab="Alleles", ylab="Profile Matches", col="green")
abline(v = 19, col = "sky blue") #indicates number of alleles where the partial match was highest
dev.off()



##### second trial #####
sample(13,7,F)
#2nd, 1st, 6th, 3rd, 13th, 9th, and 7th loci respectively
sample(2,7,T)
#1st, 2nd, 2nd, 2nd, 2nd, 2nd, and 2nd alleles respectively

#genotype controlProfile except with the corresponding loci and alleles dropped out
dropoutVec7b = c(9,0,0,25,6,0,8,10,17,15,15,0,12,0,12,11,15,0,13,12,11,12,14,18,39,0)
dropoutProfile7b = array(dropoutVec7b, dim=c(2,13))
dropoutProfile7b

#entering the SECOND dropout genotype with SEVEN alleles into function directMatch
directMatch(dropoutProfile7b, genos, ngenos)

#yielded ONE partial match at 19 out of 26 alleles matching; low chance of false positive
#yielded ONE partial match at 15 out of 26 alleles matching
#yielded TWENTY-SIX partial matches at 14 out of 26 alleles matching

sevenB = c(4,71,680,3437,11336,25234,39762,44851,37170,22698,10337,3352,877,163,26,1,0,0,0,1,0,0,0,0,0,0,0)
plot(sevenB, main = "Trial Two's Dropout at Seven Alleles", xlab="Alleles", ylab="Profile Matches", col="purple")
abline(v = 19, col = "sky blue") #indicates number of alleles where the partial match was highest

pdf(file="Trial Two's Dropout at Seven Alleles.pdf")
plot(sevenB, main = "Trial Two's Dropout at Seven Alleles", xlab="Alleles", ylab="Profile Matches", col="purple")
abline(v = 19, col = "sky blue") #indicates number of alleles where the partial match was highest
dev.off()



##### third trial #####
sample(13,7,F)
#5th, 13th, 4th, 10th, 7th, 2nd, and 9th loci respectively
sample(2,7,T)
#2nd, 2nd, 2nd, 1st, 1st, 2nd, and 2nd alleles respectively

#genotype controlProfile except with the corresponding loci and alleles dropped out
dropoutVec7c = c(9,14,23,0,6,8,8,0,17,0,15,17,0,11,12,11,15,0,0,12,11,12,14,18,39,0)
dropoutProfile7c = array(dropoutVec7c, dim=c(2,13))
dropoutProfile7c

#entering the THIRD dropout genotype with SEVEN alleles into function directMatch
directMatch(dropoutProfile7c, genos, ngenos)

#yielded ONE partial match at 19 out of 26 alleles matching; low chance of false positive
#yielded NINE partial matches at 14 out of 26 alleles matching

sevenC = c(1,106,829,3924,12723,26630,40715,44712,36336,21230,9249,2837,605,93,9,0,0,0,0,1,0,0,0,0,0,0,0)
plot(sevenC, main = "Trial Three's Dropout at Seven Alleles", xlab="Alleles", ylab="Profile Matches", col="purple")
abline(v = 19, col = "sky blue") #indicates number of alleles where the partial match was highest

pdf(file="Trial Three's Dropout at Seven Alleles.pdf")
plot(sevenC, main = "Trial Three's Dropout at Seven Alleles", xlab="Alleles", ylab="Profile Matches", col="purple")
abline(v = 19, col = "sky blue") #indicates number of alleles where the partial match was highest
dev.off()



##### fourth trial #####
sample(13,7,F)
#9th, 4th, 6th, 8th, 11th, 13th, and 7th loci respectively
sample(2,7,T)
#2nd, 1st, 2nd, 1st, 2nd, 2nd, and 1st alleles respectively

#genotype controlProfile except with the corresponding loci and alleles dropped out
dropoutVec7d = c(9,14,23,25,6,8,0,10,17,15,15,0,0,11,0,11,15,0,13,12,11,0,14,18,39,0)
dropoutProfile7d = array(dropoutVec7d, dim=c(2,13))
dropoutProfile7d

#entering the FOURTH dropout genotype with SEVEN alleles into function directMatch
directMatch(dropoutProfile7d, genos, ngenos)

#yielded ONE partial match at 16 out of 26 alleles matching; low chance of false positive
#yielded NINE partial matches at 14 out of 26 alleles matching

sevenD = c(11,248,1851,7390,19722,35054,44647,41667,28102,14076,5432,1460,278,52,9,0,0,0,0,1,0,0,0,0,0,0,0)
plot(sevenD, main = "Trial Four's Dropout at Seven Alleles", xlab="Alleles", ylab="Profile Matches", col="purple")
abline(v = 19, col = "sky blue") #indicates number of alleles where the partial match was highest

pdf(file="Trial Four's Dropout at Seven Alleles.pdf")
plot(sevenD, main = "Trial Four's Dropout at Seven Alleles", xlab="Alleles", ylab="Profile Matches", col="purple")
abline(v = 19, col = "sky blue") #indicates number of alleles where the partial match was highest
dev.off()



##### fifth trial #####
sample(13,7,F)
#13th, 10th, 5th, 9th, 12th, 3rd, and 6th loci respectively
sample(2,7,T)
#2nd, 2nd, 2nd, 2nd, 2nd, 1st, and 2nd alleles respectively

#genotype controlProfile except with the corresponding loci and alleles dropped out
dropoutVec7e = c(9,14,23,25,0,8,8,10,17,0,15,0,12,11,12,11,15,0,13,0,11,12,14,0,39,0)
dropoutProfile7e = array(dropoutVec7e, dim=c(2,13))
dropoutProfile7e

#entering the FIFTH dropout genotype with SEVEN alleles into function directMatch
directMatch(dropoutProfile7e, genos, ngenos)

#yielded ONE partial match at 16 out of 26 alleles matching; low chance of false positive
#yielded ONE partial matches at 15 out of 26 alleles matching
#yielded TWELVE partial matches at 14 out of 26 alleles matching

sevenE = c(7,164,1261,5605,16491,31275,43142,43386,31896,17293,6873,2067,473,53,12,1,0,0,0,1,0,0,0,0,0,0,0)
plot(sevenE, main = "Trial Five's Dropout at Seven Alleles", xlab="Alleles", ylab="Profile Matches", col="purple")
abline(v = 19, col = "sky blue") #indicates number of alleles where the partial match was highest

pdf(file="Trial Five's Dropout at Seven Alleles.pdf")
plot(sevenE, main = "Trial Five's Dropout at Seven Alleles", xlab="Alleles", ylab="Profile Matches", col="purple")
abline(v = 19, col = "sky blue") #indicates number of alleles where the partial match was highest
dev.off()



########## SIMULATING DROPOUT AT 10 ALLELES ##########

##### first trial #####
sample(13,10,F)
#9th, 11th, 10th, 6th, 7th, 2nd, 8th, 13th, 12th, and 3rd loci respectively
sample(2,10,T)
#2nd, 1st, 1st, 1st, 1st, 2nd, 1st, 1st, 2nd, and 1st alleles respectively

#genotype controlProfile except with the corresponding loci and alleles dropped out
dropoutVec10a = c(9,14,23,0,0,8,8,10,17,15,0,17,0,11,0,11,15,0,0,12,0,12,14,0,0,33.2)
dropoutProfile10a = array(dropoutVec10a, dim=c(2,13))
dropoutProfile10a

#entering the FIRST dropout genotype with TEN alleles into function directMatch
directMatch(dropoutProfile10a, genos, ngenos)

#yielded ONE partial match at 16 out of 26 alleles matching; low chance of false positive
#yielded FIVE partial matches at 13 out of 26 alleles matching

tenA = c(56,810,4573,14870,31276,44377,45301,32861,17159,6621,1686,358,46,5,0,0,1,0,0,0,0,0,0,0,0,0,0)
plot(tenA, main = "Trial One's Dropout at Ten Alleles", xlab="Alleles", ylab="Profile Matches", col="gray")
abline(v = 16, col = "pink") #indicates number of alleles where the partial match was highest

pdf(file="Trial One's Dropout at Ten Alleles.pdf")
plot(tenA, main = "Trial One's Dropout at Ten Alleles", xlab="Alleles", ylab="Profile Matches", col="gray")
abline(v = 16, col = "pink") #indicates number of alleles where the partial match was highest
dev.off()



##### second trial #####
sample(13,10,F)
#8th, 1st, 6th, 10th, 9th, 12th, 4th, 13th, 5th and 3rd loci respectively
sample(2,10,T)
#2nd, 2nd, 1st, 1st, 2nd, 1st, 2nd, 2nd, 2nd, and 1st alleles respectively

#genotype controlProfile except with the corresponding loci and alleles dropped out
dropoutVec10b = c(9,0,23,25,0,8,8,0,17,0,0,17,12,11,12,0,15,14,0,12,11,12,0,18,39,0)
dropoutProfile10b = array(dropoutVec10b, dim=c(2,13))
dropoutProfile10b

#entering the SECOND dropout genotype with TEN alleles into function directMatch
directMatch(dropoutProfile10b, genos, ngenos)

#yielded ONE partial match at 16 out of 26 alleles matching; low chance of false positive
#yielded TWO partial matches at 14 out of 26 alleles matching

tenB = c(14,208.1557,6756,18588,34645,45900,43198,29034,14027,4737,1130,179,24,2,0,1,0,0,0,0,0,0,0,0,0,0)
plot(tenB, main = "Trial Two's Dropout at Ten Alleles", xlab="Alleles", ylab="Profile Matches", col="gray")
abline(v = 16, col = "pink") #indicates number of alleles where the partial match was highest

pdf(file="Trial Two's Dropout at Ten Alleles.pdf")
plot(tenB, main = "Trial Two's Dropout at Ten Alleles", xlab="Alleles", ylab="Profile Matches", col="gray")
abline(v = 16, col = "pink") #indicates number of alleles where the partial match was highest
dev.off()



##### third trial #####
sample(13,10,F)
#5th, 3rd, 7th, 8th, 6th, 13th, 1st, 4th, 10th, and 2nd loci respectively
sample(2,10,T)
#1st, 1st, 2nd, 1st, 2nd, 2nd, 2nd, 1st, 1st, and 1st alleles respectively

#genotype controlProfile except with the corresponding loci and alleles dropped out
dropoutVec10c = c(9,0,0,25,0,8,0,10,0,15,15,0,12,0,0,11,15,14,0,12,11,12,14,18,39,0)
dropoutProfile10c = array(dropoutVec10c, dim=c(2,13))
dropoutProfile10c

#entering the THIRD dropout genotype with TEN alleles into function directMatch
directMatch(dropoutProfile10c, genos, ngenos)

#yielded ONE partial match at 16 out of 26 alleles matching; low chance of false positive
#yielded EIGHT partial matches at 13 out of 26 alleles matching

tenC = c(29,377,2519,9884,24128,40533,47012,38752,23137,9814,3045,649,112,8,0,0,1,0,0,0,0,0,0,0,0,0,0)
plot(tenC, main = "Trial Three's Dropout at Ten Alleles", xlab="Alleles", ylab="Profile Matches", col="gray")
abline(v = 16, col = "pink") #indicates number of alleles where the partial match was highest

pdf(file="Trial Three's Dropout at Ten Alleles.pdf")
plot(tenC, main = "Trial Three's Dropout at Ten Alleles", xlab="Alleles", ylab="Profile Matches", col="gray")
abline(v = 16, col = "pink") #indicates number of alleles where the partial match was highest
dev.off()



##### fourth trial #####
sample(13,10,F)
#2nd, 8th, 5th, 7th, 12th, 11th, 6th, 10th, 13th, and 3rd loci respectively
sample(2,10,T)
#2nd, 1st, 1st, 1st, 1st, 1st, 2nd, 2nd, 1st, and 1st alleles respectively

#genotype controlProfile except with the corresponding loci and alleles dropped out
dropoutVec10d = c(9,14,23,0,0,8,8,10,0,15,15,0,0,11,0,11,15,14,13,0,0,12,0,18,0,33.2)
dropoutProfile10d = array(dropoutVec10d, dim=c(2,13))
dropoutProfile10d

#entering the FOURTH dropout genotype with TEN alleles into function directMatch
directMatch(dropoutProfile10d, genos, ngenos)

#yielded ONE partial match at 16 out of 26 alleles matching; low chance of false positive
#yielded ONE partial matches at 14 out of 26 alleles matching
#yielded FIVE partial matches at 13 out of 26 alleles matching

tenD = c(80,1011,5285,16394,32431,44528,43754,31435,16487,6376,1804,355,53,5,1,0,1,0,0,0,0,0,0,0,0,0,0)
plot(tenD, main = "Trial Four's Dropout at Ten Alleles", xlab="Alleles", ylab="Profile Matches", col="gray")
abline(v = 16, col = "pink") #indicates number of alleles where the partial match was highest

pdf(file="Trial Four's Dropout at Ten Alleles.pdf")
plot(tenD, main = "Trial Four's Dropout at Ten Alleles", xlab="Alleles", ylab="Profile Matches", col="gray")
abline(v = 16, col = "pink") #indicates number of alleles where the partial match was highest
dev.off()



##### fifth trial #####
sample(13,10,F)
#2nd, 6th, 11th, 7th, 5th, 9th, 1st, 12th, 4th and 8th loci respectively
sample(2,10,T)
#1st, 1st, 2nd, 2nd, 1st, 2nd, 2nd, 1st, 2nd, and 1st alleles respectively

#genotype controlProfile except with the corresponding loci and alleles dropped out
dropoutVec10e = c(9,0,0,25,6,8,8,0,0,15,0,17,12,0,0,11,15,0,13,12,11,0,0,18,39.0,33.2)
dropoutProfile10e = array(dropoutVec10e, dim=c(2,13))
dropoutProfile10e

#entering the FIFTH dropout genotype with TEN alleles into function directMatch
directMatch(dropoutProfile10e, genos, ngenos)

#yielded ONE partial match at 16 out of 26 alleles matching; low chance of false positive
#yielded ONE partial matches at 14 out of 26 alleles matching
#yielded EIGHT partial matches at 13 out of 26 alleles matching

tenE = c(22,419,2906,10945,26384,41924,46875,37372,21177,8816 ,2545,542,63,8,1,0,1,0,0,0,0,0,0,0,0,0,0)
plot(tenE, main = "Trial Five's Dropout at Ten Alleles", xlab="Alleles", ylab="Profile Matches", col="gray")
abline(v = 16, col = "pink") #indicates number of alleles where the partial match was highest

pdf(file="Trial Five's Dropout at Ten Alleles.pdf")
plot(tenE, main = "Trial Five's Dropout at Ten Alleles", xlab="Alleles", ylab="Profile Matches", col="gray")
abline(v = 16, col = "pink") #indicates number of alleles where the partial match was highest
dev.off()



########## SIMULATING DROPOUT AT 13 ALLELES ##########

##### first trial #####
sample(13,13,F)
#7th, 2nd, 8th, 10th, 1st, 6th, 11th, 12th, 9th, 5th, 4th, 3rd, and 13th loci respectively
sample(2,13,T)
#1st, 2nd, 1st, 2nd, 2nd, 1st, 1st, 1st, 2nd, 2nd, 1st, 2nd, and 2nd alleles respectively

#genotype controlProfile except with the corresponding loci and alleles dropped out
dropoutVec13a = c(9,0,23,0,6,0,0,10,17,0,0,17,0,11,0,11,15,0,13,0,0,12,0,18,39,0)
dropoutProfile13a = array(dropoutVec13a, dim=c(2,13))
dropoutProfile13a

#entering the FIRST dropout genotype with THIRTEEN alleles into function directMatch
directMatch(dropoutProfile13a, genos, ngenos)

#yielded ONE partial match at both 13 out of 26 and 12 out of 26 alleles matching; low to medium chance of false positive
#yielded SEVEN partial matches at 11 out of 26 alleles matching

thirteenA = c(1226,8038,23903,42065,49168,40044,22658,9521,2713,595,60,7,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0)
plot(thirteenA, main = "Trial One's Dropout at Thirteen Alleles", xlab="Alleles", ylab="Profile Matches", col="black")
abline(v = 13, col = "orange") #indicates number of alleles where the partial match was highest

pdf(file="Trial One's Dropout at Thirteen Alleles.pdf")
plot(thirteenA, main = "Trial One's Dropout at Thirteen Alleles", xlab="Alleles", ylab="Profile Matches", col="black")
abline(v = 13, col = "orange") #indicates number of alleles where the partial match was highest
dev.off()



##### second trial #####
sample(13,13,F)
#10th, 3rd, 6th, 4th, 1st, 5th, 9th, 13th, 7th, 11th, 2nd, 12th, and 8th loci respectively
sample(2,13,T)
#2nd, 2nd, 1st, 1st, 1st, 2nd, 1st, 2nd, 2nd, 2nd, 2nd, 2nd, and 2nd alleles respectively

#genotype controlProfile except with the corresponding loci and alleles dropped out
dropoutVec13b = c(0,14,23,0,6,0,0,10,17,0,0,17,12,0,12,0,0,14,13,0,11,0,14,0,39,0)
dropoutProfile13b = array(dropoutVec13b, dim=c(2,13))
dropoutProfile13b

#entering the SECOND dropout genotype with THIRTEEN alleles into function directMatch
directMatch(dropoutProfile13b, genos, ngenos)

#yielded ONE partial match at 13 out of 26 alleles matching; low to medium chance of false positive
#yielded EIGHT partial matches at 11 out of 26 alleles matching

thirteenB = c(743,5782,19973,40374,51316,43656,25236,9820,2569,480,42,8,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0)
plot(thirteenB, main = "Trial Two's Dropout at Thirteen Alleles", xlab="Alleles", ylab="Profile Matches", col="black")
abline(v = 13, col = "orange") #indicates number of alleles where the partial match was highest

pdf(file="Trial Two's Dropout at Thirteen Alleles.pdf")
plot(thirteenB, main = "Trial Two's Dropout at Thirteen Alleles", xlab="Alleles", ylab="Profile Matches", col="black")
abline(v = 13, col = "orange") #indicates number of alleles where the partial match was highest
dev.off()



##### third trial #####
sample(13,13,F)
#3rd, 5th, 7th, 9th, 10th, 4th, 1st, 8th, 6th, 2nd, 12th, 11th, and 13th loci respectively
sample(2,13,T)
#1st, 1st, 1st, 2nd, 1st, 1st, 2nd, 1st, 2nd, 2nd, 2nd, 2nd, and 2nd alleles respectively

#genotype controlProfile except with the corresponding loci and alleles dropped out
dropoutVec13c = c(9,0,23,0,0,8,0,10,0,15,15,0,0,11,0,11,15,0,0,12,12,0,18,0,39,0)
dropoutProfile13c = array(dropoutVec13c, dim=c(2,13))
dropoutProfile13c

#entering the THIRD dropout genotype with THIRTEEN alleles into function directMatch
directMatch(dropoutProfile13c, genos, ngenos)

#yielded ONE partial match at both 13 out of 26 and 12 out of 26 alleles matching; low to medium chance of false positive
#yielded EIGHTEEN partial matches at 11 out of 26 alleles matching

thirteenC = c(321,3252,13186,30759,46761,47494,33660,16938,6025,1346,238,18,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0)
plot(thirteenC, main = "Trial Three's Dropout at Thirteen Alleles", xlab="Alleles", ylab="Profile Matches", col="black")
abline(v = 13, col = "orange") #indicates number of alleles where the partial match was highest

pdf(file="Trial Three's Dropout at Thirteen Alleles.pdf")
plot(thirteenC, main = "Trial Three's Dropout at Thirteen Alleles", xlab="Alleles", ylab="Profile Matches", col="black")
abline(v = 13, col = "orange") #indicates number of alleles where the partial match was highest
dev.off()



##### fourth trial #####
sample(13,13,F)
#9th, 12th, 13th, 3rd, 4th, 10th, 11th, 5th, 2nd, 6th, 7th, 1st, and 8th loci respectively
sample(2,13,T)
#1st, 2nd, 1st, 2nd, 2nd, 2nd, 1st, 2nd, 2nd, 2nd, 1st, 2nd, and 2nd alleles respectively

#genotype controlProfile except with the corresponding loci and alleles dropped out
dropoutVec13d = c(9,0,23,0,6,0,8,0,17,0,15,0,0,11,12,0,0,14,13,0,0,12,14,0,0,33.2)
dropoutProfile13d = array(dropoutVec13d, dim=c(2,13))
dropoutProfile13d

#entering the FOURTH dropout genotype with THIRTEEN alleles into function directMatch
directMatch(dropoutProfile13d, genos, ngenos)

#yielded ONE partial match at both 13 out of 26 and 12 out of 26 alleles matching; low to medium chance of false positive
#yielded FIFTEEN partial matches at 11 out of 26 alleles matching

thirteenD = c(407,3612,14276,32294,47483,47027,32377,15638,5433,1242,194,15,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0)
plot(thirteenD, main = "Trial Four's Dropout at Thirteen Alleles", xlab="Alleles", ylab="Profile Matches", col="black")
abline(v = 13, col = "orange") #indicates number of alleles where the partial match was highest

pdf(file="Trial Four's Dropout at Thirteen Alleles.pdf")
plot(thirteenD, main = "Trial Four's Dropout at Thirteen Alleles", xlab="Alleles", ylab="Profile Matches", col="black")
abline(v = 13, col = "orange") #indicates number of alleles where the partial match was highest
dev.off()



##### fifth trial #####
sample(13,13,F)
#11th, 8th, 5th, 9th, 6th, 2nd, 3rd, 1st, 4th, 10th, 7th, 12th, and 13th loci respectively
sample(2,13,T)
#2nd, 1st, 1st, 2nd, 1st, 2nd, 1st, 1st, 1st, 1st, 1st, 1st, and 1st alleles respectively

#genotype controlProfile except with the corresponding loci and alleles dropped out
dropoutVec13e = c(0,14,23,0,0,8,0,10,0,15,0,17,0,11,0,11,15,0,0,12,11,0,0,18,0,33.2)
dropoutProfile13e = array(dropoutVec13e, dim=c(2,13))
dropoutProfile13e

#entering the FIFTH dropout genotype with THIRTEEN alleles into function directMatch
directMatch(dropoutProfile13e, genos, ngenos)

#yielded ONE partial match at 13 out of 26 alleles matching; low to medium chance of false positive
#yielded TWENTY partial matches at 11 out of 26 alleles matching

thirteenE = c(310,3146,12625,30126,46288,47551,34314,17593,6249,1556,221,20,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0)
plot(thirteenE, main = "Trial Five's Dropout at Thirteen Alleles", xlab="Alleles", ylab="Profile Matches", col="black")
abline(v = 13, col = "orange") #indicates number of alleles where the partial match was highest

pdf(file="Trial Five's Dropout at Thirteen Alleles.pdf")
plot(thirteenE, main = "Trial Five's Dropout at Thirteen Alleles", xlab="Alleles", ylab="Profile Matches", col="black")
abline(v = 13, col = "orange") #indicates number of alleles where the partial match was highest



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


