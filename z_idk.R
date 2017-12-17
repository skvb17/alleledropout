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
controlProfile <- genos[30732,,]




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




#function for one allele dropout at a random loci SOFTCODE
  #input: selected genotype (profile)
  #purpose: changes matrix so that it changes chosen allele to 0 effectively dropping an allele out
  #output: genotype with one random allele dropped out

controlVec = c(9,14,23,25,6,8,8,10,17,15,15,17,12,11,12,11,15,14,13,12,11,12,14,18,39.0,33.2)
controlProfile = array(controlVec, dim=c(2,13))
controlProfile

dropout1 = function(profile) {
  locidrop <- sample(13,1,F)
  alledrop <- sample(2,1,T) 
  profile[alledrop, locidrop] <- 0

return(profile)
}

dropout1(controlProfile) #works!!


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
controlProf05a


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


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


##### JACKY'S CODE FOR DOING ALLELE DROPOUTS AT 1 TO 13 ALLELES ####
CONTROL <- controlProfile
for (j in 1:13){ 
  dropout7 = function(DDprofile) {
    LOCIDROP <- sample(13,j,F)
    ALLEDROP <- sample(2,j,T) 
    for (i in 1:j)
      if(DDprofile[ALLEDROP[i], LOCIDROP[i]] <- 0) {
        print(DDprofile)
      }
    return(DDprofile)
  }
  print(data.frame(dropout7(CONTROL)))
  #  print(directMatch(CONTROL, genos, ngenos)) 
}

#make a for loop for directmatch later

#trying to assign arrays to variables???? need helpppp!
for (j in 1:13) { 
  dropout = function(DDprofile) {
    LOCIDROP <- sample(13,j,F)
    ALLEDROP <- sample(2,j,T) 
    for (i in 1:j)
      if(DDprofile[ALLEDROP[i], LOCIDROP[i]] <- 0) {
        print(DDprofile)
      }
    return(DDprofile)
  }
  
  print(dropout(CONTROL))
  #print(directMatch(CONTROL, genos, ngenos)) 
}


#for (k in 1:) {
directMatch = function(array, database, nProfiles) {
  directmatches = rep(0,27)
  for (i in 1:nProfiles) {
    directmatches0 = alleleMatch(array,database[i,,])
    directmatches[directmatches0 + 1] = directmatches[directmatches0 + 1] + 1
  }
  
  return(directmatches)
}


#STUFF WE NEED TO FIGURE OUT
  #how to use our information to input it into directMatch
  


result = NULL
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
  print(directMatch(result, genos, ngenos))
}








