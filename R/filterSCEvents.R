# filterSCEvents
#
# Custom function extracted from CustomFunctions not Annotated New 250131.Rmd
# (lines 45-55). This is also exported by the CustomFunctionsSept25 R package;
# included here as a sourced script so the analysis .Rmd is self-contained.

filterSCEvents <- function(sce, n_cells = NULL)
  
{
  
test = table(sce$sample_id)
test = test[test>n_cells]
test = as.data.frame(test)
sce = filterSCE(sce, sample_id %in% test$Var1)
  
return(sce)
}
