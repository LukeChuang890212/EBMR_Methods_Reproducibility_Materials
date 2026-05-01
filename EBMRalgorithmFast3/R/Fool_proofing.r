check_data = function(y_names, data) {
  if(sum(is.na(data[y_names])) > 0){
    data$r = ifelse(is.na(data[y_names]), 0, 1)
  }else{
    if(!("r" %in% colnames(data))){
      stop("there must be a column of response indicator with column name \"r\".")
    }
  }
  return(data)
}
