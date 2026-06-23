# create_superscript_column
# 
# Convert merging labels to superscript-aware character column (lines 1280-1323 of source .Rmd).
# Migrated from CMV CyTOF Figures David.Rmd as part of repository reorganisation
# (see CMV_paper_analysis.Rmd / CMV_extra_analyses.Rmd / CMV_code_quarantine.Rmd).

create_superscript_column <- function(sce, cell_id_column, search_terms = NULL, replacement_terms = NULL) {

  # Convert factor to character if necessary
  cell_ids <- as.character(sce[[cell_id_column]])
  
  # Check for any NA values and warn if found
  if (any(is.na(cell_ids))) {
    warning("There are NA values in the cell IDs. They will not be modified.")
  }

  print("Original Cell IDs:")
  print(unique(cell_ids))  # Debug: Check the original unique values

  # Convert minus and plus to superscript versions by default
  cell_ids <- gsub("-", "^'−'*", cell_ids)   # Superscript minus sign from hyphen
  cell_ids <- gsub("\\+", "^'+'*", cell_ids) # Superscript plus sign
  cell_ids <- gsub(" ", "~", cell_ids)       # Convert spaces to tildes

  print("Cell IDs after default superscripting:")
  print(unique(cell_ids))  # Debug: Check after minus and plus transformation

  # Apply additional search and replace terms if provided
  if (!is.null(search_terms)) {
    for (i in seq_along(search_terms)) {
      print(paste0("Superscripting term: ", search_terms[[i]], " to: ^'", replacement_terms[[i]], "'*"))
      cell_ids <- gsub(search_terms[[i]], paste0("^'", replacement_terms[[i]], "'*"), cell_ids)
    }
  }

  # Remove any asterisks (*) at the end of the string
  cell_ids <- gsub("\\*+$", "", cell_ids)

  print("Final Cell IDs:")
  print(unique(cell_ids))  # Debug: Check final transformed IDs

  # Add the new superscripted column to the sce object
  sce[[paste0(cell_id_column, "_superscript")]] <- cell_ids

  return(sce)
}


