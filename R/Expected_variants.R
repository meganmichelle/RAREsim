#' RAREsim
#'
#' Simulate rare variant genetic data
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
  
  ### columns need to be names lower and upper
  #### want  an error if  ONE isn't true!
  if((colnames(mac_bin_prop)[1] == 'Lower') == FALSE | (colnames(mac_bin_prop)[2]  == 'Upper') == FALSE
     | (colnames(mac_bin_prop)[3]  == 'Prop') ==  FALSE){
    stop('mac_bin_prop needs to have column names Lower, Upper, and  Prop, respectively')
  }
  
  ## mac needs  to have numeric values
  if((is.numeric(mac_bin_prop$Lower) ==  FALSE) | (is.numeric(mac_bin_prop$Upper) == FALSE)
     | (is.numeric(mac_bin_prop$Prop) == FALSE)){
    stop('The columns of mac_bin_prop need to be numberic')
  }

  mac_bin_prop$Expected_var <- Total_num_var*mac_bin_prop$Prop
  
  mac <- mac_bin_prop[,c(1,2,4)]
  
  return(mac)
}

