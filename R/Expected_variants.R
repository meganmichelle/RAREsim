#' Expected Variants
#'
#' This function combines the Number of Variants and AFS functions to produces the expected number of variants per Kb in each MAC bin
#' The output from the afs and nvariants functions can be used here
#' 
#' 
#' @param mac_bin_prop The MAC bins to use, with lower and upper boundaries defined. Only define for rare variants
#' 
#' @param Total_num_var estimated  number of variants  in the  region
#'
#' @return data frame with the MAC bins and expected variants
#'
#' @author Megan M Null, \email{megan.null@ucdenver.edu}
#' 
#' @keywords RAREsim
#'
#' @export


Expected_variants <- function(Total_num_var,mac_bin_prop){
  

  if(is.numeric(Total_num_var) == FALSE){
    stop('Error: Total_num_var needs to be numeric')
  }
  
  if((colnames(mac_bin_prop)[1] == 'Lower') == FALSE | (colnames(mac_bin_prop)[2]  == 'Upper') == FALSE
     | (colnames(mac_bin_prop)[3]  == 'Prop') ==  FALSE){
    stop('mac_bin_prop needs to have column names Lower, Upper, and  Prop, respectively')
  }
  
  ## mac needs  to have numeric values in each column
  if(is.numeric(mac_bin_prop$Prop) == FALSE){
    stop('The column of Proportion of variants is required to be numeric')
  }

  # Multiply the proportion in each bin by the total number of variants
  mac_bin_prop$Expected_var <- Total_num_var*mac_bin_prop$Prop
  
  # Subset to the columns we want (removing the original proporiton)
  mac <- mac_bin_prop[,c(1,2,4)]
  
  # Return the new data frame with the expected number of variants per bin
  return(mac)
}

