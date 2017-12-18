#getting database "genos"
getwd()
setwd("~/Desktop/Shannel-project")
load("~/Desktop/Shannel-project/Rori's Database.Rdata")


#getting allele frequency table
freq = read.table("AFs_AfA.dat.txt") #reads in data
t(freq)



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

#getting random profile for experiment and turning it into an array
dim(genos)  # 200K profiles 2 allelles  13 loci
dim(genos)[1] # 200k
ngenos = dim(genos)[1]
sample(ngenos,1) -> PROFILE  # choosing # of RANDOM profile(s)
PROFILE
genos[PROFILE,,] -> CONTROL # this require a manual insert of the choosen profile # 
CONTROL

#entering the control genotype into function directMatch
directMatch(CONTROL, genos, ngenos) 


# Creating a empty DataFrame name DROP13
DROP13 <- matrix(0, ncol = 27, nrow = 13)
DROP13 <- data.frame(DROP13)
rownames(DROP13) <- c("Drop1","Drop2", "Drop3", "Drop4","Drop5", "Drop6", "Drop7","Drop8", "Drop9", "Drop10","Drop11","Drop12", "Drop13")

# set Variables to NULL for the forLoops
result =NULL
directmatchResult=NULL

############# 
#####
# Result for ONE profile with 1 to 13 alleles drop

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
                  #print(unlist(directMatch(result, genos, ngenos)))  #remove # to have the 
                  }

# Creating a empty DataFrame name DROP_DIFF
DROP_DIFF <- matrix(0, ncol = 27, nrow = 13)
DROP_DIFF <- data.frame(DROP_DIFF)
rownames(DROP_DIFF) <- c("Drop1","Drop2", "Drop3", "Drop4","Drop5", "Drop6", "Drop7","Drop8", "Drop9", "Drop10","Drop11","Drop12", "Drop13")


############# 
#####
# Result for DIFFERENT profiles with 1 to 13 alleles drop

result =NULL
directmatchResult=NULL
  
  for (j in c(1:13)){ 
    
    sample(ngenos,1) -> PROFILE  # choosing # of RANDOM profile(s)
    PROFILE
    genos[PROFILE,,] -> CONTROL # this require a manual insert of the choosen profile # 
    CONTROL
    
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
    #print(unlist(directMatch(result, genos, ngenos)))
    directmatchResult = unlist(directMatch(result, genos, ngenos)) -> DROP_DIFF[j,]
    print(unlist(directMatch(result, genos, ngenos)))
  }

