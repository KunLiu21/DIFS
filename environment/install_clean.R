# Fresh Linux installation, with no author library or package cache.
if (as.character(getRversion())!='4.4.1') stop('Use R 4.4.1 for this recipe')
lib <- Sys.getenv('R_LIBS_USER')
if (!nzchar(lib)) stop('Set R_LIBS_USER to a new, empty directory')
dir.create(lib,recursive=TRUE,showWarnings=FALSE)
if (length(list.files(lib,all.files=TRUE,no..=TRUE))) stop('R_LIBS_USER must initially be empty')
.libPaths(c(normalizePath(lib),.Library),include.site=FALSE)
if (length(.libPaths())!=2L || .libPaths()[1]!=normalizePath(lib)) stop('Unexpected library search path')
dir.create('clean-room/evidence',recursive=TRUE,showWarnings=FALSE)
writeLines(capture.output(list(R=R.version.string,libPaths=.libPaths(),
  initial_packages=installed.packages()[,c('Package','Version','LibPath','Priority')],
  empty_user_library=TRUE)),'clean-room/evidence/before-install.txt')
cran <- 'https://packagemanager.posit.co/cran/__linux__/jammy/2024-09-30'
options(repos=c(CRAN=cran),timeout=1200,Ncpus=2,
  BioC_mirror='https://bioconductor.posit.co')
install.packages(c('BiocManager','remotes'),lib=lib)
BiocManager::install(version='3.19',lib=lib,ask=FALSE,update=FALSE)
# First obtain CRAN dependencies from one dated snapshot, then select the
# historical top-level versions from CRAN's version archive without upgrading.
install.packages(c('Seurat','diptest','matrixStats','igraph','magrittr'),lib=lib)
pins <- c(Matrix='1.7-0',SeuratObject='5.0.2',Seurat='4.4.0',
          diptest='0.77-2',matrixStats='1.4.1',igraph='2.0.3')
for (p in names(pins)) {
  actual <- tryCatch(as.character(utils::packageVersion(p)),error=function(e) '')
  if (actual!=pins[[p]]) remotes::install_version(p,version=pins[[p]],lib=lib,
    repos='https://cloud.r-project.org',dependencies=FALSE,upgrade='never')
}
options(repos=c(CRAN=cran))
BiocManager::install(c('SC3','FEAST','DuoClustering2018','ExperimentHub'),
  version='3.19',lib=lib,ask=FALSE,update=FALSE,Ncpus=2)
for (p in names(pins)) if (as.character(utils::packageVersion(p))!=pins[[p]])
  stop('Pinned version mismatch: ',p)
for (p in c('SC3','FEAST','DuoClustering2018','ExperimentHub','Seurat','diptest'))
  if (!requireNamespace(p,quietly=TRUE)) stop('Required package did not load: ',p)
if (as.character(packageVersion('SC3'))!='1.32.0' ||
    as.character(packageVersion('FEAST'))!='1.12.0') stop('Bioconductor target versions differ')
ip <- as.data.frame(installed.packages(fields='Repository')[,c('Package','Version','LibPath','Repository')],stringsAsFactors=FALSE)
write.csv(ip,'clean-room/evidence/installed-packages.csv',row.names=FALSE)
writeLines(capture.output(sessionInfo()),'clean-room/evidence/install-session.txt')
cat('PASS: fresh package installation completed.\n')
