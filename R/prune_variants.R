#' Theoretically prune the variants
#'
#' This function theoretically prunes the variants to match what is expected
#' The current allele counts from the haplotype file also needs to be provided
#' The output will be used to prune the haplotype files, as seen in the example on github.
#'
#' @param MAC vector of allele count per variant that corresponds with the haplotype file. 
#' Line k of the MAC file needs to correspond with the variant on line k in the haplotype/legend files
#'
#' @param expected dataframe of expected number of variants per bin - as well as bin boundaries
#'
#' @return A list with two data frames: ToRemove (variants that will have all alternate alleles pruned) and ToChange (variants that will have a subset of alternate alleles removed), Both have 3 columns: line number, current MAC, new MAC
#'
#'
#' @export
#' 
#' @importFrom stats runif


prune_variants <- function(MAC, expected){
  
  # check for correct column names
  if( ((colnames(expected)[1] == 'Lower') & (colnames(expected)[2] == 'Upper') 
     & (colnames(expected)[3] == 'Expected_var') ) == FALSE){
    stop('Error: Column names from expected should be "Lower", "Upper", and "Expected_var"')
  }
  
  # check that expected is  all numeric
  if((is.numeric(expected$Lower) ==  FALSE) | (is.numeric(expected$Upper) == FALSE)
     | (is.numeric(expected$Expected_var) == FALSE)){
    stop('The columns of expected need to be numeric')
     }

  # MAC bins need to be ordered from smallest to largest
  if(is.unsorted(expected$Lower)==TRUE){
  stop('The MAC bins need to be ordered from smallest to largest')
  }
  
  # MAC bins need to be non-overlapping and exhaustive
  temp <- c()
  for(i in 1:nrow(expected)){
    tmp <- c(expected$Lower[i]: expected$Upper[i])
    temp <- c(temp, tmp)
  }
  check <- c(expected$Lower[1]:expected$Upper[nrow(expected)])
  if(identical(temp, check) == FALSE){
    print('Warning: The MAC bins should be disjoint and exhaustive with respect to rare allele counts.
          Check your MAC bins.')
  }
  
  rem_all<-c() # store variants that will have all alternate alleles removed
  change_all <-  c()  # store variants that will have a subset of alternate alleles removed
  
  # check that there are more simulated rare alleles than expected
  if(length(which(MAC$V1>0))<sum(expected$Expected_var)){
    print('Error, not enough rare alleles')
    print('Check that >2000 individuals are being simulated and the size of the region is correct')
  }
  
  MAC$num <-  1:nrow(MAC) # Identify the line each variant is on
  colnames(MAC)[1] <- 'V1' # consistent columnn names for the function
  
  for(k in nrow(expected):1){ # loop through the bins, largest to smallest
    tmp <-  MAC[which(MAC$V1 >= expected$Lower[k] & MAC$V1 <= expected$Upper[k]),] # pull the variants with the MAC bin of interest
    
    tmp <- as.data.frame(tmp)
    
    colnames(tmp)[1] <- 'V1'
    
    n1 <- nrow(tmp) # number of simulated variants in the MAC bin
    
    expect <- expected$Expected_var[k] # number of expected  variants in the MAC bin

    if(n1>0){ # when there are simulated variants in the bin
      p <- expect/n1 # calculate expected over simulated
    }else{ # when there are not any simulated variants in the bin
      if(expect < 3){next} # skip if we expect less than 3
      if(expect >=3){ # if we expect more than 3 variants in the bin
        p <- expect # p is just the number we expect
      }
    }
    
    # When we have less expected variants than simulated in the MAC bin
    if(p <= 1){ 
      tmp$rd <- runif(n1) # take a random draw from the uniform(0,1) distribution for each variant
      rem <- tmp[c(which(tmp$rd >= p)),] # remove the variants >= p (so keep variants with probability p)
      rem_all <- rbind(rem_all, rem) # create a list of all variants to remove
    }
    
    # When we have more expected variants than simulated in the MAC bin
    if(p>1){
        if((expect-n1) <3){ # skip the bin if the number of simulated variants is within 3 of expected
        next}
      n2 <- expect  - n1 ### n2 is how many more we need
      p1 <- n2/nrow(rem_all)  ### this is the proportion from rem_all to instead only remove a subset of variants
      
      if(p1 > 1){ # if there are not enough pruned variants to add, give a warning and add all
        print(paste0('Warning: The Bin with MAC ', expected$Lower[k], ':', expected$Upper[k],
                     ' has less than the expected number of minor alleles'))
        to_change <- rem_all # add all variants that were previously pruned
        mac_range <-  c(expected$Lower[k]:expected$Upper[k]) # list the possible new MACs
        to_change$new_MAC <- sample(mac_range, nrow(to_change),
                                    replace = TRUE) # randomly select a new MAC
        rem_all <- c() # all the pruned variants are now removing a subset of alternate alleles
        change_all<- rbind(change_all, to_change) # bind the list of variants to the full list
        next
      }
      
      rem_all$rd <- runif(nrow(rem_all)) # new random draw for all of the variants to prune
      to_change <- rem_all[which(rem_all$rd  <= p1),] # select variants to remove a subset with probability p1
      mac_range <-  c(expected$Lower[k]:expected$Upper[k]) # list the possible new MACs
      to_change$new_MAC <- sample(mac_range, nrow(to_change),
                                  replace = TRUE) # randomly select a new MAC
      
      rem_all <- rem_all[which(rem_all$rd > p1),] # remove the variants that will have a subset of alternated alleles returned to reference from the list of completely pruned variants 
      change_all<- rbind(change_all, to_change) # create a complete list of all the variants to remove a subset of alternate alleles
    }
    
  } # end of the loop through the bins
  
  rem_all <- rem_all[,c(2,1)] # remove the random draw and reorder columns
  colnames(rem_all) <- c('line', 'Current_mac') # rename colmns to be informative
  rem_all$New_mac <- 0 # All removed variants will have a new MAC = 0
  
  change_all <- change_all[,c(2,1,4)] # reorder the columns to make rem_all
  if(is.null(change_all)==FALSE){ # If there are variants that need a subset of alternate alleles returned to reference
  colnames(change_all) <- colnames(rem_all) # rename the columns to make rem_all
  }
  
  return(list(ToRemove = rem_all, ToChange = change_all))
}






