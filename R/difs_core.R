# Compatibility entry point for stage I only; full DIFS is in R/difs.R.
.core_file <- sys.frame(1)$ofile
source(file.path(dirname(.core_file), "difs.R"))
.core_engine <- difs_load_engine(full=FALSE)
for (.name in c("difs_gate_threshold", "difs_gate_genes", "difs_stage1_ranking",
                "dip_pvalue_lookup_vec", "difs_null_tables"))
  assign(.name, .core_engine[[.name]], envir=globalenv())
rm(.core_file, .name)
