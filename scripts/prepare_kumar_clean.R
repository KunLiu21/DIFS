# Run the documented data entry point using new, explicitly empty hub caches.
source('R/difs.R')
difs_require(c('ExperimentHub','AnnotationHub'))
cache <- normalizePath('clean-room',mustWork=TRUE)
eh <- file.path(cache,'ExperimentHub'); ah <- file.path(cache,'AnnotationHub')
for (d in c(eh,ah)) {
  dir.create(d,recursive=TRUE,showWarnings=FALSE)
  if (length(list.files(d,all.files=TRUE,no..=TRUE))) stop('Hub cache is not empty: ',d)
}
if (file.exists('example/input/Kumar.rds')) stop('Prepared input already exists')
Sys.setenv(EXPERIMENT_HUB_CACHE=eh,ANNOTATION_HUB_CACHE=ah)
ExperimentHub::setExperimentHubOption('CACHE',eh)
ExperimentHub::setExperimentHubOption('ASK',FALSE)
AnnotationHub::setAnnotationHubOption('CACHE',ah)
AnnotationHub::setAnnotationHubOption('ASK',FALSE)
writeLines(c('prepared_input_present_before=FALSE','hub_caches_empty_before=TRUE',
  paste0('ExperimentHub_version=',packageVersion('ExperimentHub')),
  paste0('DuoClustering2018_version=',packageVersion('DuoClustering2018')),
  paste0('ExperimentHub_cache=',eh),paste0('AnnotationHub_cache=',ah)),
  'clean-room/evidence/fresh-data-check.txt')
source('example/prepare_kumar.R')
cat('PASS: public Kumar download and documented preparation completed.\n')
