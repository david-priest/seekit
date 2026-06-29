# sceget.R — CATALYST-free rewrite promoted from dev/catalyst_quarantine/.
# Body verbatim from the project-vendored copy; only the CATALYST namespace
# shims (CATALYST:::.* internals, bare accessors, the asNamespace('CATALYST')
# hack) were rewritten to the package's .wl_* internals (R/wl_internals.R).
# 2026-06 seekit migration of the CMV CyTOF pipeline.
# sceget
#
# Custom function extracted from CustomFunctions not Annotated New 250131.Rmd
# (lines 710-747). This is also exported by the CustomFunctionsSept25 R package;
# included here as a sourced script so the analysis .Rmd is self-contained.

sceget <- function(sce, k = "merging1", meta = c("patient_id"), meta2 = "sample_id ~ cluster_id",ncells = NULL, merging_col = NULL)
  
{
  
  # Filter the sce by number of cells if desired
if(!is.null(ncells)){
  sce <- filterSCEvents(sce, n_cells = ncells)
}
  
if (!is.null(merging_col)) { 
cluster_ids <- sce[[k]] 
} else { 
k <- .wl_check_k(sce, k) 
cluster_ids <- .wl_cluster_ids(sce, k) 
}
  
# Make a data frame containing cluster abundances
ns <- table(cluster_id = cluster_ids, sample_id = .wl_sample_ids(sce))
fq <- prop.table(ns, 2) * 100
df <- as.data.frame(fq)
m <- match(df$sample_id, sce$sample_id)

# Extract any other metadata into the df data frame (replace the variables in c() with your own choices).
for (i in meta) df[[i]] <- sce[[i]][m]

#dfout <<- df

# dcast to get rows as patient ids
df = dcast(df, meta2, value.var = "Freq")

#dfout2 <<- df
# put it in the same order as ei(sce)  !! Important for diffcyt!  I'm not doing this for now!
#sampids <- ei(sce)

#df <- df[match(sampids$sample_id, df$sample_id),]

return(df)
}
