# Compare newly computed results with the archived executed example.
source('R/difs.R')
original <- 'example/validated/kumar-20260925'
fresh <- 'example/output-clean'
old <- read.csv(file.path(original,'cell_labels.csv'),colClasses='character')
new <- read.csv(file.path(fresh,'cell_labels.csv'),colClasses='character')
stopifnot(setequal(old$cell,new$cell))
new <- new[match(old$cell,new$cell),]
stopifnot(identical(old$truth,new$truth))
old_genes <- read.csv(file.path(original,'selected_genes.csv'))$gene
new_genes <- read.csv(file.path(fresh,'selected_genes.csv'))$gene
comparison <- data.frame(
  cells=nrow(new),gene_sets_identical=setequal(old_genes,new_genes),
  gene_order_identical=identical(old_genes,new_genes),
  shared_genes=length(intersect(old_genes,new_genes)),
  final_partition_ARI=unname(difs_metrics(setNames(old$predicted,old$cell),
    setNames(new$predicted,new$cell))['ARI']),
  preliminary_partition_ARI=unname(difs_metrics(setNames(old$preliminary,old$cell),
    setNames(new$preliminary,new$cell))['ARI']))
write.csv(comparison,'clean-room/evidence/old-versus-fresh.csv',row.names=FALSE)
print(comparison)
if (!comparison$gene_sets_identical || comparison$final_partition_ARI < 1-1e-12 ||
    comparison$preliminary_partition_ARI < 1-1e-12)
  stop('Fresh run differs from the archived selections/partitions; inspect before claiming full agreement')
cat('PASS: fresh data/install run reproduces the archived genes and partitions.\n')
