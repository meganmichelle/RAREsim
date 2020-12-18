#' RAREsim
#'
#' Simulate rare variant genetic data
#'
#' @param con1 original name of the file
#'
#' @param con2 output file
#'
#' @param repl_lines replacement lines
#'
#' @param line_nos line numbers to be  replaces
#'
#' @return Text stating what has been accomlished
#'
#' @author Megan M Null, \email{megan.null@ucdenver.edu}
#' 
#' @keywords RAREsim
#'
#' @export
#' 

update_hap_files<-function(con1, con2, repl_lines, line_nos){
  f1 = file(con1, open="r")
  f2 = file(con2, open="w")
  l = readLines(f1, warn=FALSE)
  for (i  in 1:length(l)) {
    if (i %in% line_nos) {
      ind = match(i, line_nos)
      write(repl_lines[ind], f2)
    } else {
      write(l[i], f2)
    }
  }
  
  print(paste("Processed", length(l), "records"))
  print(paste("Replaced", length(line_nos), "lines in original file", con1))
  print(paste("New file is", con2))
  
  close(f1)
  close(f2)
}

