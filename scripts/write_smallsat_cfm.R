#--------------------------------------------------------------------------------------------
# Write out Confusion matrices to CSV
#--------------------------------------------------------------------------------------------
write_smallsat_cfm = function(out.file) {
#' Confusion Matrices to CSV
#'
#' @param out.file `string`; name and location (path) of the output file
#'
#' @return `.csv` file exported to the `out.file` location 
#' @export
#'
#' @examples

  
  # Error matrix -- MLHC (not maj)
  PA = error.mtx.mlhc$byClass[,3] # Pos Pred Value >> not sure if this is correct?
  UA = error.mtx.mlhc$byClass[,4] # Neg Pred Value >> not sure if this is correct?
  mat = error.mtx.mlhc$table
  order.txt = c("gt","rt","yt","hb","bg","sh")
  order = match(order.txt,colnames(mat))
  mat = mat[c(order),c(order)]
  mat = rbind(mat,colSums(mat),UA[order[1:6]])
  mat = cbind(mat,rowSums(mat),c(PA[order[1:6]],NA,NA))
  rownames(mat) = c(order.txt,"sum","UA")
  colnames(mat) = c(order.txt,"sum","PA")
  mat
  info.mat = paste("error.mtx.mlhc$overall"," Overall Acc:",round(as.numeric(error.mtx.mlhc$overall[1]*100),3),"%")
  write.table(info.mat, file.path(exportsPATH, out.file), append = FALSE, sep = ',', col.names = F, row.names = F)
  write.table(mat, file.path(exportsPATH, out.file ), append = TRUE, sep = ',', col.names = T, row.names = TRUE)
  write.table(c("-------"), file.path(exportsPATH, out.file), append = TRUE, sep = ',', col.names = TRUE, row.names = TRUE)
  #- - - - - - - - - - - - - - - - - 
  # Ollofson Error matrix -- MLHC (not maj)
  mat = olof.table.mlhc$matrix
  order.txt = c("gt","rt","yt","hb","bg","sh","sum")
  order = match(order.txt,colnames(mat))
  mat = mat[c(order),c(order)]
  mat = rbind(mat,c(olof.table.mlhc$PA[order[1:6]],NA))
  mat = cbind(mat,c(olof.table.mlhc$UA[order[1:6]],NA,NA))
  rownames(mat) = c(order.txt,"UA")
  colnames(mat) = c(order.txt,"PA")
  mat
  info.mat = paste("olof.table.mlhc$matrix, Overall Acc:",round(olof.table.mlhc$OA*100,3),"%")
  write.table(info.mat, file.path(exportsPATH,  out.file), append = TRUE, sep = ',', col.names = F, row.names = F)
  write.table(mat, file.path(exportsPATH,  out.file), append = TRUE, sep = ',', col.names = c(order.txt,"PA"), row.names = TRUE)
  write.table(c("-------"), file.path(exportsPATH,  out.file), append = TRUE, sep = ',', col.names = TRUE, row.names = TRUE)
  #--------------------------------------------------------------------------------------------
  # Error matrix -- rf 
  PA = error.mtx.rf$byClass[,3] 
  UA = error.mtx.rf$byClass[,4] 
  mat = error.mtx.rf$table
  order.txt = c("gt","rt","yt","hb","bg","sh")
  order = match(order.txt,colnames(mat))
  mat = mat[c(order),c(order)]
  mat = rbind(mat,colSums(mat),UA[order[1:6]])
  mat = cbind(mat,rowSums(mat),c(PA[order[1:6]],NA,NA))
  rownames(mat) = c(order.txt,"sum","UA")
  colnames(mat) = c(order.txt,"sum","PA")
  mat
  info.mat = paste("error.mtx.rf$overall"," Overall Acc:",round(as.numeric(error.mtx.rf$overall[1]*100),3),"%")
  write.table(info.mat, file.path(exportsPATH,  out.file), append = T, sep = ',', col.names = F, row.names = F)
  write.table(mat, file.path(exportsPATH, out.file), append = TRUE, sep = ',', col.names = T, row.names = TRUE)
  write.table(c("-------"), file.path(exportsPATH,  out.file), append = TRUE, sep = ',', col.names = TRUE, row.names = TRUE)
  #- - - - - - - - - - - - - - - - - 
  # Ollofson Error matrix -- rf
  mat = olof.table.rf$matrix
  order.txt = c("gt","rt","yt","hb","bg","sh","sum")
  order = match(order.txt,colnames(mat))
  mat = mat[c(order),c(order)]
  mat = rbind(mat,c(olof.table.rf$PA[order[1:6]],NA))
  mat = cbind(mat,c(olof.table.rf$UA[order[1:6]],NA,NA))
  rownames(mat) = c(order.txt,"UA")
  colnames(mat) = c(order.txt,"PA")
  mat
  info.mat = paste("olof.table.rf$matrix, Overall Acc:",round(olof.table.rf$OA*100,3),"%")
  write.table(info.mat, file.path(exportsPATH,  out.file), append = TRUE, sep = ',', col.names = F, row.names = F)
  write.table(mat, file.path(exportsPATH,  out.file), append = TRUE, sep = ',', col.names = c(order.txt,"PA"), row.names = TRUE)
  write.table(c("-------"), file.path(exportsPATH,  out.file), append = TRUE, sep = ',', col.names = TRUE, row.names = TRUE)
  #--------------------------------------------------------------------------------------------
  # Error matrix -- svm
  PA = error.mtx.svm$byClass[,3] 
  UA = error.mtx.svm$byClass[,4]
  mat = error.mtx.svm$table
  order.txt = c("gt","rt","yt","hb","bg","sh")
  order = match(order.txt,colnames(mat))
  mat = mat[c(order),c(order)]
  mat = rbind(mat,colSums(mat),UA[order[1:6]])
  mat = cbind(mat,rowSums(mat),c(PA[order[1:6]],NA,NA))
  rownames(mat) = c(order.txt,"sum","UA")
  colnames(mat) = c(order.txt,"sum","PA")
  mat
  info.mat = paste("error.mtx.svm$overall"," Overall Acc:",round(as.numeric(error.mtx.svm$overall[1]*100),3),"%")
  write.table(info.mat, file.path(exportsPATH,  out.file), append = T, sep = ',', col.names = F, row.names = F)
  write.table(mat, file.path(exportsPATH,  out.file), append = TRUE, sep = ',', col.names = T, row.names = TRUE)
  write.table(c("-------"), file.path(exportsPATH,  out.file), append = TRUE, sep = ',', col.names = TRUE, row.names = TRUE)
  #- - - - - - - - - - - - - - - - - 
  # Ollofson Error matrix -- svm 
  mat = olof.table.svm$matrix
  order.txt = c("gt","rt","yt","hb","bg","sh","sum")
  order = match(order.txt,colnames(mat))
  mat = mat[c(order),c(order)]
  mat = rbind(mat,c(olof.table.svm$PA[order[1:6]],NA))
  mat = cbind(mat,c(olof.table.svm$UA[order[1:6]],NA,NA))
  rownames(mat) = c(order.txt,"UA")
  colnames(mat) = c(order.txt,"PA")
  mat
  info.mat = paste("olof.table.svm$matrix, Overall Acc:",round(olof.table.svm$OA*100,3),"%")
  write.table(info.mat, file.path(exportsPATH,  out.file), append = TRUE, sep = ',', col.names = F, row.names = F)
  write.table(mat, file.path(exportsPATH,  out.file), append = TRUE, sep = ',', col.names = c(order.txt,"PA"), row.names = TRUE)
  write.table(c("-------"), file.path(exportsPATH,  out.file), append = TRUE, sep = ',', col.names = TRUE, row.names = TRUE)
  #--WRITE VARIABLE SELECTION----------------------
  write.table(paste("vARIABLE SELECTION (mir):",varsel.treshold,"Tresholds"), file.path(exportsPATH,  out.file), append = TRUE, sep = ',', col.names = TRUE, row.names = TRUE)
  write.table(print(sort(vardat,decreasing = T)), file.path(exportsPATH,  out.file), append = TRUE, sep = ',', col.names = TRUE, row.names = TRUE)
  write.table(c("-------"), file.path(exportsPATH,  out.file), append = TRUE, sep = ',', col.names = TRUE, row.names = TRUE)
  
  #---------------------------------------------------------------------------------------------------
}
