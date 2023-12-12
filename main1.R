

library(predictR)
outputDir <- '/home/pgsalome/R/toolbx/organized_folders/NTCP_TCP/output_NEW'
features_folder <- '/home/pgsalome/R/toolbx/organized_folders/NTCP_TCP/clidmeraddos_features/ricci-dose/lem1'
outcomecsv <- paste('/home/pgsalome/R/toolbx/organized_folders/ricci_outcomes.csv',sep='')
# datacsv <- "/home/pgsalome/R/toolbx/organized_folders/NTCP_TCP/clidmeraddos_features/ricci/features_cli_ohe.csv"
# fsSurv(datacsv,outcomecsv,"os",outputDir,100)
# List of outcome types
outcome_types <- c("os", "pfs", "tox")
iter = 250
run_parallel = TRUE
# Get list of feature files
feature_files <- list.files(features_folder, pattern = "features_.*\\.csv$", full.names = TRUE)


# Main loop
for (file in feature_files) {
  parts <- unlist(strsplit(get_mdroiname(file), "_"))
  mod <- parts[1]
  roi <- parts[2]

  for (outcometype in outcome_types) {
    outputDir_md <- paste(outputDir, "/feature_selection_results2/", outcometype, "/", mod, "_", roi, sep = '')

    # Create output directory if it doesn't exist
    if (!dir.exists(outputDir_md)) {
      dir.create(outputDir_md, recursive = TRUE)
    }

    # Call fsSurv function
    fsSurv(file, outcomecsv, outcometype, outputDir_md,iter,run_parallel) # Assuming '2' is a fixed parameter in your context
  }
}

# fsSignature(features_folder, outputDir, freq = 500)
# fsImpute(outputDir)

#
#
# cohort1 = 'ricci_firstorder'
# cohort = 'ricci'
# result_folder = paste('/home/pgsalome/R/toolbx/organized_folders/NTCP_TCP/featuresignificance_results/',cohort,sep='')
# features_folder = paste('/home/pgsalome/R/toolbx/organized_folders/NTCP_TCP/clidmeraddos_features/',cohort,sep='')
# outputDir = paste('/home/pgsalome/R/toolbx/organized_folders/NTCP_TCP/modroi_signatures2/',cohort,sep='')
#
#
# outcomecsv = paste('/home/pgsalome/R/toolbx/organized_folders/',cohort,'_outcomes.csv',sep='')
# toolbxdir = '/home/pgsalome/R/toolbx'
# outcometype = 'tox'
# name = 'letdds'
# impute  = FALSE
# top3 = FALSE
# freq = 900
# p = 0.05
# sls = 0.05
#
# dir.create(paste(outputDir,'/',outcometype,'/Imputed',sep=''),recursive = TRUE)
# dir.create(paste(outputDir,'/',outcometype,'/Complete',sep=''),recursive = TRUE)
#
#
# #### leave impute at false because you want to compare uch vs imp
# fsSignature(features_folder,toolbxdir,outcomecsv,outcometype,outputDir,to_include,name,sls,p,freq,top3,impute = impute)
