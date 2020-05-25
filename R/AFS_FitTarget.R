
Fit_AFS <- function(prob, N, per_rv){ ### only works with N>2200
  ### prob is the probability of each bin (9 probabilities)
  ### N is the number of samples
  n05 <- floor(N*2*.005) ### MAF = 0.5%
  n1 <- floor(N*2*.01) ### MAF = 1%
  
  
  
  mac <- as.data.frame(matrix(nrow = 7, ncol = 2))
  colnames(mac) <- c('Start', 'End')
  mac$Start <- c(1:3,6,11,21,(n05+1))
  mac$End <- c(1:2,5,10,20,n05,n1)
  
  mac$prob <- (as.numeric(as.character(prob))) #### this will need to be cleaned up so it  always works
  names(mac)[3] <- 'prob'
  #print(mac)
  
  
  c1 <- 1:n1 #### MACs to use in the function
  
  #suppressMessages()
  calc_prob_LS<-function(tune){
    ### calculate b
    alpha_mac <- 1/((c1+tune[2])^(tune[1]))
    b <- per_rv/sum(alpha_mac) #### now we have b!
    b_mac <- b*alpha_mac
    all <- 0
    for(i in 1:nrow(mac)){
      E <- sum(b_mac[mac$Start[i]:mac$End[i]])
      O <- mac$prob[i]
      c <- (E-O)^2
      all <- all+c
    }
    return(all)
  }  
  
  tune <- c(1,0) #### start with the function 1/x
  S <- slsqp(tune, fn = calc_prob_LS, hin = hin.tune )
  b <- per_rv/sum(1/((c1+S$par[2])^(S$par[1])))
  
  return(c(S$par[2], S$par[1], b))
  #print(c(S$par[2], S$par[1], b))
  
}
