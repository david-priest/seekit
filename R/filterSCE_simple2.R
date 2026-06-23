# filterSCE_simple2
#
# Custom function extracted from CustomFunctions not Annotated New 250131.Rmd
# (lines 792-808). This is also exported by the CustomFunctionsSept25 R package;
# included here as a sourced script so the analysis .Rmd is self-contained.

filterSCE_simple2 <- function(sce, filterColumn, filterValues, exclude = FALSE) {
  # Get the condition
  if(exclude) {
    condition <- !sce[[filterColumn]] %in% filterValues
  } else {
    condition <- sce[[filterColumn]] %in% filterValues
  }

  # Filter the sce object
  sce_filtered <- sce[, condition]
  
  # Drop unused levels in all factor columns in the colData
  is_factor <- sapply(colData(sce_filtered), is.factor)
  colData(sce_filtered)[is_factor] <- lapply(colData(sce_filtered)[is_factor], droplevels)

  return(sce_filtered)
}
