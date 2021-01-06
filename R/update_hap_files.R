#' Edit the (upzipped) haplotype files
#' 
#' After creating new pruned variants with the *create_new_mac* function,
#' edit the unzipped haplotypes. A new haplotype file is created.
#' This needs to be done before variants are removed to ensure the line numbers align correctly. 
#'
#'
#' @param con1 name of the original unzipped haplotype file
#'
#' @param con2 name of the new haplotype file, also unzipped
#'
#' @param repl_lines replacement lines - obtained from readLines of the file to replace
#'
#' @param line_nos line numbers to be replaced
#'
#' @return R output is text stating what has been accomplished. The new files have been written.
#'
#'
#' @export
#' 

update_hap_files<-function(con1, con2, repl_lines, line_nos){
  
  f1 = file(con1, open="r") # the original haplotype file
  f2 = file(con2, open="w") # the name of the new haplotype file
  l = readLines(f1, warn=FALSE)
  for (i  in 1:length(l)) { # loop through the lines of the file
    if (i %in% line_nos) {
      ind = match(i, line_nos)
      write(repl_lines[ind], f2)  # replace the new variants
    } else {
      write(l[i], f2) # do nothing if the variant doesn't need updated
    }
  }
  
  print(paste("Processed", length(l), "records"))
  print(paste("Replaced", length(line_nos), "lines in original file", con1))
  print(paste("New file is", con2))
  
  close(f1)
  close(f2)
}

