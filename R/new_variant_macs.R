#' Edits the haplotype file for variants that require a subset of alternate alleles returned to reference
#'
#' The *new_variant_macs* function edits the variants that require a subset of alternate alleles returned to reference.
#' The variants that will be edited should first extracted from the haplotype file and the new MAC for each variant determined.
#' See the example scripts for how to extract the variants, and the *pruning_info* function to determine which
#' variants should have a subset of alternate alleles returned to reference.
#'
#' @param original_variants A data frame with the haplotypes of the variants the require a subset of the alternate alleles to be returned to reference
#'
#' @param new_MAC The new MAC for each variant. The variants counts need to match the line number of original_haps.
#'
#' @return New dataframe with the edited variants to insert into the complete haplotype files
#'
#'
#' @export
#' 

new_variant_macs<-function(original_variants, new_MAC){
  ### create data frame with all reference alleles to edit with new alternate alleles:
  updated_mac  <-  data.frame(matrix(0, ncol = ncol(original_variants),
                                     nrow = nrow(original_variants))) 
  
  for(i in 1:nrow(original_variants)){ ## loop through each variant
    
    k <- new_MAC[i] ## extract the number of alternate alleles to keep
    aa <- which(original_variants[i,]>0)  ### create of list of the haps (columns) with MAC >0
    aa_n <- sample(aa, k) #### sample k alternate alleles that will remain alternate alleles
    
    updated_mac[i,aa_n] <- 1 ### change the selected haplotypes to alternate alleles
    
    print(i)
  }
  
  return(updated_mac)
}

