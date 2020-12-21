#' RAREsim
#'
#' Edit the necessary variants
#'
#' @param original_haps
#'
#' @param new_MAC
#'
#' @return New dataframe with the edited variants
#'
#' @author Megan M Null, \email{megan.null@ucdenver.edu}
#' 
#' @keywords RAREsim
#'
#' @export
#' 

new_variant_macs<-function(original_variants, new_MAC){
  new_variants <- original_variants
  for(i in 1:nrow(original_variants)){
    temp <- original_variants[i,]
    k <- new_MAC[i]
    ### create of list of the haps (columns) with MAC >0
    aa <- which(temp>0)
    #### sample k alternate alleles that will remain alternate alleles
    aa_n <- sample(aa, k)
    updated_mac  <-  rep(0,length(temp)) 
    updated_mac[aa_n] <- 1
    
    ### use the new row
    new_variants[i,] <- updated_mac
    print(i)
  }
  
  return(new_variants)
}
  

