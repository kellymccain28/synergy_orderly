# devtools::install_github('mrc-ide/hipercow@mrc-6733')
# devtools::install_github('mrc-ide/hipercow@mrc-6733', subdir = "drivers/dide")
# install.packages(
#   c("hipercow.dide"),
#   repos = c("https://mrc-ide.r-universe.dev",
#             "https://cloud.r-project.org"))
library(hipercow)

cyphr::data_admin_init(".")


x <- runif(10)
key <- cyphr::data_key()
cyphr::encrypt(saveRDS(x, "data.rds"), key)

hipercow_provision(method = 'pkgdepends', refs = c('cyphr', 'mrc-ide/hipercow@mrc-6733'))

# generate a keypair 
dide_generate_keypair(update = TRUE)

t <- task_create_expr(cyphr::data_request_access())
task_wait(t)

# request access
cyphr::data_admin_list_requests() #don't see anything from this access request
cyphr::data_admin_authorise(yes = TRUE)

task_status(t)
task_log_show(t)

t <- task_create_expr(
  cyphr::decrypt(readRDS("data.rds"), cyphr::data_key()))
task_result(t)
