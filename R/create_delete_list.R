#' Create file for effiecient pruning of variants
#' 
#' @details 
#' Given a list of variants, create_delete_list will write a file that specified which row/variants to delete from the haplotype and legend files.
#' The written file will be executed with a sed or gsed command in bash.
#'
#' @param variants_to_remove list of the line  number of variants to remove
#' 
#' @param name name of the file that will be written. Default to 'List2delete.sed'
#'
#' @return R output is text stating what has been accomlished. A file has been written.
#' 
#' @examples 
#' create_delete_list(1:10)
#' #In bash: sed -i -f List2delete.sed [file to delete]
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
  
  output.file = file(paste0(name,'.sed'), 'wb')
  # write the code, in a file ending in .sed
  write.table(to_remove, file = output.file, row.names = FALSE,
              col.names = FALSE, quote = FALSE)
  close(output.file)
  
  print(paste0('File written called ', name, '.sed'))
  print(paste0('To implement, run: sed -i -f ',
               name, '.sed [file to delete]'))
  
}

