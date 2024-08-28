#--------------------------------------------------------------------------------------------
# Write out metrics to CSV
#--------------------------------------------------------------------------------------------
write_smallsat_metrics = function(out.file) {
#' Metrics to CSV
#'
#' @param out.file `string`; name and location (path) of the output file
#'
#' @return `.csv` file file exported to the `out.file` location 
#' @export
#'
#' @examples
  
  header = c("img","Total_area(km2)","Total_NA","Total_mask","Total_data",
             "OA_mlhc","OA_mlhc_olof","OA_rf","OA_rf_olof","OA_svm","OA_svm_olof",
             names(n.class.rf.arr),names(n.class.rf.arr),names(n.class.rf.arr))
  OA_metrics = c(img_id,area.img,area.na,area.mask,area.data,
                 as.numeric(error.mtx.mlhc$overall[1]),olof.table.mlhc$OA,
                 as.numeric(error.mtx.rf$overall[1]),olof.table.rf$OA,
                 as.numeric(error.mtx.svm$overall[1]),olof.table.svm$OA,
                 as.numeric(c(n.class.mlhc.arr[],n.class.rf.arr[],n.class.svm.arr[])))
  write.table(rbind(header,OA_metrics), file.path(exportsPATH, out.file), append = FALSE, sep = ',', col.names = F, row.names = F)
}