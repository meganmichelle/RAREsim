#' RAREsim
#'
#' Simulate rare variant genetic data
#'
#' @param variants_to_remove list of the line  number of variants to remove
#' 
#' @param name name of the file that will be written
#'
#' @return Text stating what has been accomlished
#'
#' @author Megan M Null, \email{megan.null@ucdenver.edu}
#' 
#' @keywords RAREsim
#'
#' @export
#' 

create_delete_list<-function(variants_to_remove, name=NULL){
  
  if(is.null(name)){
    name = 'List2delete'
  }

  to_remove <- paste(variants_to_remove, collapse = 'd; ')
  to_remove <- paste0(to_remove,'d')

  write.table(to_remove, file = paste0(name,'.sed'), row.names = FALSE, 
              col.names = FALSE, quote = FALSE)
  
  print(paste('File written in working directory called', name))
  
}

