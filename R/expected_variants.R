#' Combines Number of Variants and AFS functions
#'
#' This function combines the Number of Variants and AFS functions to produces the expected number of variants per Kb in each MAC bin
#'
#'
#' @param mac_bin_prop The MAC bins to use, with three columns: Lower, Upper, and Prop.
#' Lower and Upper define the MAC bins boundaries and Prop is the proportion of variants in each respective bin.
#'  Only define for rare variants
#'
#' @param Total_num_var estimated total number of variants in the region of interest
#'
#' @return data frame with the MAC bins and expected variants
#'
#' @examples
#'  data('afs_afr')
#'  mac <- afs_afr[,c(1:2)]
#'  expected_variants(Total_num_var = 19.029*nvariant(pop='AFR', N = 8128),
#'  mac_bin_prop = afs(mac_bins = mac, pop = 'AFR'))
#'
#' @export


expected_variants <- function(Total_num_var,mac_bin_prop){
    
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

