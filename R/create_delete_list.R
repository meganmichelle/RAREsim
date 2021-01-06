#' Create file for effiecient pruning of variants
#'
#' @param variants_to_remove list of the line  number of variants to remove
#' 
#' @param name name of the file that will be written. Default to 'List2delete.sed'
#'
#' @return R output is text stating what has been accomlished. A file has been written.
#' 
#'
#' @export
#' 
#' @importFrom utils write.table

create_delete_list<-function(variants_to_remove, name=NULL){
  
  if(is.null(name)){
    name = 'List2delete'
  }

  # bind all the numbers together in a way that can be read by the sed command
  to_remove <- paste(variants_to_remove, collapse = 'd; ')
  
  # add another 'd' at the end for the command
  to_remove <- paste0(to_remove,'d')

  # write the code, in a file ending in .sed
  write.table(to_remove, file = paste0(name,'.sed'), row.names = FALSE, 
              col.names = FALSE, quote = FALSE)
  
  print(paste('File written called', name))
  print(paste0('To implement, run: sed -i -f ',
  name, '.sed [New_haps.haps]'))
  
}

