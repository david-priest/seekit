# sce_meta
#
# Custom function extracted from CustomFunctions not Annotated New 250131.Rmd
# (lines 1327-1346). This is also exported by the CustomFunctionsSept25 R package;
# included here as a sourced script so the analysis .Rmd is self-contained.

sce_meta <- function(sce, metafields,md_u)
  
{

  #https://stackoverflow.com/questions/35636315/replace-values-in-a-dataframe-based-on-lookup-table
  df = as.data.frame(sce@colData$sample_id)
  lookup = md_u
  
  for (i in 1:length(metafields)){
  
  field <- metafields[i]  
  #print(paste("Now updating: ",field))
  new <- df
  new[] <- lapply(df, function(x) lookup[[field]][match(x, lookup$sample_id)])
  sce[[field]] = as.factor(new$`sce@colData$sample_id`)
    
  }

  return(sce)
}
