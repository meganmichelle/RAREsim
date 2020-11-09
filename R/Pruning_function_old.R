#' RAREsim
#'
#' Simulate rare variant genetic data
#'
#' @param hap_file_name Name of the haplotype file to prune (gzipped: end in .gz)
#'
#' @param MAC vector of allele count per variant that corresponds with the haplotype file
#'
#' @param expected dataframe of expected number of variants per bin - as well as bin boundaries
#'
#' @param sed_or_gsed  Specify if gsed is loaded onto computer instead of sed
#'
#' @return Number of variants pruned
#'
#' @author Megan M Null, \email{megan.null@ucdenver.edu}
#' 
#' @keywords RAREsim
#'
#' @export
#' 
#' @importFrom utils write.table
#' @importFrom stats runif

Pruning_function <- function(hap_file_name, MAC, expected, sed_or_gsed = NULL){
  
  rem_all<-c()
  change_all <-  c() 
  if(length(which(MAC$V1>0))<sum(expected$Expected_var)){
    print('Error, not enough rare alleles')
  }
  MAC$num <-  1:nrow(MAC)
  for(k in nrow(expected):1){
      tmp <-  MAC[which(MAC$V1 >= expected$Lower[k] & MAC$V1 <= expected$Upper[k]),]
      tmp <- as.data.frame(tmp)
      colnames(tmp)[1] <- 'V1'
        ##### if there are not any SNVs in the bin, tmp is all NA's...
      n1 <- nrow(tmp)
      expect <- expected$Expected_var[k]
        #print(expect)
      if(n1>0){
      p <- expect/n1
      }else{
          if(expect < 1){next}
          if(expect >=1){
            p <- expect
          }
        }
        if(p <= 1){
          
          rd <- runif(n1)
          rem <- tmp[c(which(rd >= p)),]
          rem_all <- rbind(rem_all, rem)
        }
        if(p>1){
          if(k==7){
            if((expect-n1) >2){
              print(c('Error!',n1, expect, (n1-expect)))} ## instead just skip
            next}
          n2 <- expect  - n1 ### n2 is how many more we need
          p1 <- n2/nrow(rem_all)  ### this is the proportion from rem_all
          rem_all$rd <- runif(nrow(rem_all))
          to_change <- rem_all[which(rem_all$rd  <= p1),]
          mac_range <-  c(expected$Lower[k]:expected$Upper[k])
          to_change$new_MAC <- sample(mac_range, nrow(to_change),
                                      replace = TRUE)
          ### change them within the MAC file.
          ####  this needs to change for the haplotype file!!!
          rem_all <- rem_all[which(rem_all$rd > p1),]
          change_all<- rbind(change_all, to_change)
        }
  }
  ##### here we are done with the theoretical pruning
  ##### now change the haplotypes
  ### First, the variants that change alternate alleles
  command_unzip <- paste0('gunzip ', hap_file_name)
  system(command_unzip)
  hap_file_name1 <- as.character(hap_file_name)
  hap_file_name1 <- substr(hap_file_name1, 1, (nchar(hap_file_name1)-3))
     if(length(change_all)>0){
        ##### change all: MAC, row #, draw, new MAC
        for(ch in 1:nrow(change_all)){
          
          n <- change_all$num[ch] #### this is the line of interest
          k <- change_all$new_MAC[ch]
        
        get_line <- paste0('sed -n \'', n, 'p;',
                      (n+1),'q\' ', hap_file_name1)
        ri <- system(get_line, intern = TRUE)
        re <- unlist(strsplit(ri, split=" "))
        re <- as.numeric(as.character(re))
        
        #### find which have alternate alleles
        aa <- which(re>0)
        
        #### sample
        aa_n <- sample(aa, k)
        
        #### and now make a new row to insert into the file
        ref <- seq(0, length(re)) ### not sure this works
        ref[aa_n] <- 1
        
        ##### format to write into the file
        ref1 <- paste(ref, collapse = ' ')
        
        #### write to a new file
        write.table(ref1, 'Text2add.txt', row.names = FALSE, 
                    col.names = FALSE, quote = FALSE)
        
        
        torun <- paste0('printf "%s\n" "', n, 'r Text2add.txt" w | ed -s ', 
                        hap_file_name1)
        system(torun)
        #### and remove the 'orginal' line n
        if(sed_or_gsed == 'gsed'){
        system(paste0('gsed -i \'', n, 'd\' ', hap_file_name1))
        }else{
          system(paste0('sed -i \'', n, 'd\' ', hap_file_name1))
        }
        }
     }
       
       ##### this needs to work witth the zipped file!
       todel <- paste(rem_all$num, collapse = 'd; ')
       todel <- paste0(todel,'d')
       write.table(todel, 'List2delete.sed', row.names = FALSE, 
                   col.names = FALSE, quote = FALSE)
       #### delete the entire list
       
       if(sed_or_gsed == 'gsed'){
         torun <- paste0('gsed -i -f List2delete.sed ', hap_file_name1)
       }else{
         torun <- paste0('sed -i -f List2delete.sed ', hap_file_name1)
       }
       
       system(torun)
       system(paste0('gzip ', hap_file_name1))
       return(nrow(rem_all))
}


Pruning_function2_linux <- function(hap_file_name, MAC, expected){
  if(is.character(hap_file_name)==FALSE){
    stop('hap_file_names needs to be  a character')
  }
  
  if(is.character(hap_file_name)==FALSE){
    stop('hap_file_names needs to be  a character')
  }
  if((substr(hap_file_name, (nchar(hap_file_name)-2), nchar(hap_file_name)) == '.gz') == FALSE){
    stop('A gzipped haplotype file is required. String must end in .gz')
  }
  
  if((ncol(expected)  == 3) ==FALSE){
    stop('The expected data frame requries 3 columns.')
  }
  
  if(is.numeric(MAC[,1])==FALSE){
    stop('The MAC file is required to be numeric')
  }
  
  if(is.null(sed_or_gsed)){
    sed_or_gsed <- 'sed'
  }
  
  
  rem_all<-c()
  change_all <-  c() 
  if(length(which(MAC$V1>0))<sum(expected$Expected_var)){
    print('Error, not enough rare alleles')
  }
  MAC$num <-  1:nrow(MAC)
  colnames(MAC)[1] <- 'V1'
  for(k in nrow(expected):1){
    tmp <-  MAC[which(MAC$V1 >= expected$Lower[k] & MAC$V1 <= expected$Upper[k]),]
    if(nrow(tmp) == 0){
      next
    }
    tmp <- as.data.frame(tmp)
    colnames(tmp)[1] <- 'V1'
    ##### if there are not any SNVs in the bin, tmp is all NA's...
    n1 <- nrow(tmp)
    expect <- expected$Expected_var[k]
    #print(expect)
    if(n1>0){
      p <- expect/n1
    }else{
      if(expect < 1){next}
      if(expect >=1){
        p <- expect
      }
    }
    if(p <= 1){
      tmp$rd <- runif(n1)
      rem <- tmp[c(which(tmp$rd >= p)),]
      rem_all <- rbind(rem_all, rem)
    }
    if(p>1){
      if(k==7){
        if((expect-n1) >2){
          print(c('Error!',n1, expect, (n1-expect)))} ## instead just skip
        next}
      n2 <- expect  - n1 ### n2 is how many more we need
      p1 <- n2/nrow(rem_all)  ### this is the proportion from rem_all
      rem_all$rd <- runif(nrow(rem_all))
      to_change <- rem_all[which(rem_all$rd  <= p1),]
      mac_range <-  c(expected$Lower[k]:expected$Upper[k])
      to_change$new_MAC <- sample(mac_range, nrow(to_change),
                                  replace = TRUE)
      ### change them within the MAC file.
      ####  this needs to change for the haplotype file!!!
      rem_all <- rem_all[which(rem_all$rd > p1),]
      change_all<- rbind(change_all, to_change)
    }
  }
  ##### here we are done with the theoretical pruning
  ##### now change the haplotypes
  ### First, the variants that change alternate alleles
  command_unzip <- paste0('gunzip ', hap_file_name)
  system(command_unzip)
  hap_file_name1 <- as.character(hap_file_name)
  hap_file_name1 <- substr(hap_file_name1, 1, (nchar(hap_file_name1)-3))
  
  con1 <- hap_file_name1
  con2 <- paste0(hap_file_name1, '_new.haps')
  lin_nos <- change_all$num
  
  
  if(length(change_all)>0){
    
    new <- as.data.frame(matrix(nrow=nrow(change_all),  ncol = 113770))
    ##### change all: MAC, row #, draw, new MAC
    for(ch in 1:nrow(change_all)){
      
      n <- change_all$num[ch] #### this is the line of interest
      k <- change_all$new_MAC[ch]
      
      get_line <- paste0('sed -n \'', n, 'p;',
                         (n+1),'q\' ', hap_file_name1)
      ri <- system(get_line, intern = TRUE)
      re <- unlist(strsplit(ri, split=" "))
      re <- as.numeric(as.character(re))
      
      #### find which have alternate alleles
      aa <- which(re>0)
      
      #### sample
      aa_n <- sample(aa, k)
      
      #### and now make a new row to insert into the file
      ref <- rep(0,length(re)) 
      ref[aa_n] <- 1
      
      
      #### write to a new file
      
      new[ch,] <- ref
      
      #### and remove the 'orginal' line n
    }
    
    
    f3 = file("test_full_text2add1.txt", "w")
    for (i in 1:nrow(new)){
      ref1 <- paste(new[i,], collapse = ' ')
      write(ref1, f3)
    }
    close(f3)
    
    repl_lines = readLines("test_full_text2add1.txt")
    
    system.time(m_update(con1, con2, repl_lines, line_nos  =  lin_nos))
    
  }else{
    system(paste0('cp ', con1, ' ', con2))
  }
  #### copy the original ones if nothing to change??
  
  ####  then  go in and remove the other  lines
  
  todel <- paste(rem_all$num, collapse = 'd; ')
  todel <- paste0(todel,'d')
  write.table(todel, 'List2delete.sed', row.names = FALSE, 
              col.names = FALSE, quote = FALSE)
  #### delete the entire list
  
  if(sed_or_gsed == 'gsed'){
    torun <- paste0('gsed -i -f List2delete.sed ', hap_file_name1)
  }else{
    torun <- paste0('sed -i -f List2delete.sed ', hap_file_name1)
  }
  
  system(torun)
  system(paste0('gzip ', con2))
  return(nrow(rem_all))
}


