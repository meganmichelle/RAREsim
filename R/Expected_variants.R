#' RAREsim
#'
#' Simulate rare variant genetic data
#'
#' @param temp A data frame with the first column sample size and the second variants per Kb
#'
#' @return data frame with  the bins, proportion, and expected variants
#'
#' @author Megan M Null, \email{megan.null@ucdenver.edu}
#' @seealso \code{\link{AFS_FitTarget}}
#' @keywords RAREsim
#'
#'
#' @export
#'
#'

Expected_variants <- function(alpha, beta, b, phi, omega, Ntar, Size){

  
  m <- floor(Ntar*2*0.01) ### RV MAC in the target data
  fit <- as.data.frame(matrix(nrow=1, ncol = m))
  
  for(i in 1:m){
    fit[,i] <- b/((beta+i)^alpha)
  }
  
  #### add the loop matrix here?? that might be easiest
  loop <- as.data.frame(matrix(nrow=7,  ncol = 3))
  colnames(loop) <- c('Start', 'End', 'Prop')
  loop$Start <- c(1,2,3,6,11,21,m/2)
  loop$End <- c(1,2,5,10,20,((m/2)+1),m)
  
  for(k in 1:nrow(loop)){
    loop$Prop[k] <- sum(fit[1,c(loop$Start[k]:loop$End[k])])
  }
  
  loop$Expected_var <- (Size/1000)*(phi*Ntar^omega)*loop$Prop
  
  return(loop)
}

