#' RAREsim
#'
#' Simulate rare variant genetic data
#'
#' @param prop_df data frame with 3 columns, Start, End (of MAC bin)  and proportion of variants in that MAC bin
#'
#' @param N Number of individuals in the target data
#'
#' @param p_rv percent of rare variants - just the sum of the proportions
#'
#' @return Vector of parameters - alpha, beta, and b as well as  fitted values
#'
#' @author Megan M Null, \email{megan.null@ucdenver.edu}
#' 
#' @keywords RAREsim
#'
#'
#' @export
#' @importFrom nloptr slsqp
#'

### prop is the proportion of each bin (7 proportions)
### N is the number of samples
# n05 <- floor(N*2*.005) ### MAF = 0.5%
# n1 <- floor(N*2*.01) ### MAF = 1%
# 
# 
# 
# mac <- as.data.frame(matrix(nrow = 7, ncol = 2))
# colnames(mac) <- c('Start', 'End')
# mac$Start <- c(1:3,6,11,21,(n05+1))
# mac$End <- c(1:2,5,10,20,n05,n1)

Fit_AFS <- function(prop_df, N, p_rv){ ### only works with N>2200
  
  prop_df$prop <- (as.numeric(as.character(prop_df$prop))) 
  
  hin.tune <- function(x) {
    h <- numeric(1)
    h[1] <- x[1] 
    return(h)
  }

  c1 <- 1:prop_df$End[nrow(prop_df)] #### MACs to use in the function
  
  #suppressMessages()
  calc_prob_LS<-function(tune){
    ### calculate b
    alpha_mac <- 1/((c1+tune[2])^(tune[1]))
    b <- p_rv/sum(alpha_mac) #### now we have b!
    b_mac <- b*alpha_mac
    all <- 0
    for(i in 1:nrow(prop_df)){
      E <- sum(b_mac[prop_df$Start[i]:prop_df$End[i]])
      O <- prop_df$prop[i]
      c <- (E-O)^2
      all <- all+c
    }
    return(all)
  }  
  
  tune <- c(1,0) #### start with the function 1/x
  S <- suppressMessages(slsqp(tune, fn = calc_prob_LS, hin = hin.tune ))
  b <- p_rv/sum(1/((c1+S$par[2])^(S$par[1])))
  
  ### now get the proportion for each bin?

  re <- prop_df
  re$prop <- '.'
  colnames(re)[3] <- 'Fitted_prop'
  
  fit <- as.data.frame(matrix(nrow =  1,  ncol = prop_df$End[nrow(prop_df)]))
  
  for(i in 1:prop_df$End[nrow(prop_df)]){
    fit[,i] <- b/((S$par[2]+i)^S$par[1])
  }
  
  for(i in 1:nrow(re)){
  re$Fitted_prop[i] <- sum(fit[,c(re$Start[i]:re$End[i])])
  }
  
  return(list(parameters = c(S$par[2], S$par[1], b), Fitted_results = re))
  
}
