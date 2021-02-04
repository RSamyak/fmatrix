library(fmatrix)
library(phylodyn)
library("phangorn")
library("mvtnorm")
library("ape")

root <- rprojroot::has_file("DESCRIPTION")
setwd(root$find_file("vignettes/beast_compare/"))

source("simulation_dna_for_Beast.R")

dataname <- "temp4"

true_tree <- readRDS(paste0(dataname, ".RDS"))
true_Fmat <- gen_Fmat(true_tree)


mcc_tree <- ape::read.nexus(paste0(dataname, ".mcc"))
mcc_Fmat <- gen_Fmat(mcc_tree)
is.Fmat(mcc_Fmat)

plot.new()
plotF(mcc_Fmat)



schedule <- "lin"
load(paste(dataname, schedule, "schedule.RData", sep = "_"))
lin_min.config <- track.config[which.min(track.energy.config), ]

schedule <- "log"
load(paste(dataname, schedule, "schedule.RData", sep = "_"))
log_min.config <- track.config[which.min(track.energy.config), ]

schedule <- "exp"
load(paste(dataname, schedule, "schedule.RData", sep = "_"))
exp_min.config <- track.config[which.min(track.energy.config), ]


par(mfrow = c(2,3))

exp_min.config %>%
  tree_from_matching %>% gen_Fmat %>% plotF
lin_min.config %>%
  tree_from_matching %>% gen_Fmat %>% plotF
log_min.config %>%
  tree_from_matching %>% gen_Fmat %>% plotF

energy(matching(mcc_tree))
energy(lin_min.config)
energy(log_min.config)
energy(exp_min.config)

mcc_Fmat %>% plotF

true_Fmat %>% plotF

dev.off()

## Question 0:  Is True Tree meaningful?
## Question 1:  Is the True Frechet Mean of the Posterior distribn close to the True Tree?
## Question 2:  Is our approximation of Frechet Mean close to True Frechet Mean?

## Question A:  Point estimate vs Credible sets
distance_Fmat(true_Fmat, mcc_Fmat)
distance_Fmat(true_Fmat, exp_min.config %>%
                tree_from_matching %>% gen_Fmat)
distance_Fmat(true_Fmat, log_min.config %>%
                tree_from_matching %>% gen_Fmat)
distance_Fmat(true_Fmat, lin_min.config %>%
                tree_from_matching %>% gen_Fmat)





par(mfrow = c(1,3))
out <- lapply(Fmat_list, distance_Fmat, Fmat2 = gen_Fmat(treeP))
hist(unlist(out))
out2 <- lapply(Fmat_list, distance_Fmat, Fmat2 = exp_min.config %>%
                 tree_from_matching %>% gen_Fmat)
hist(unlist(out2))
out_mcc <- lapply(Fmat_list, distance_Fmat, Fmat2 = gen_Fmat(mcc_tree))
hist(unlist(out_mcc))

subtr.list <- ape::subtrees(treeP, wait=TRUE)
child.ntip <- data.frame(node.id=(1:Nnode(treeP))+n.tip,
                         n.tip=unlist(lapply(subtr.list, Ntip)))

subtr.list_mcc <- ape::subtrees(mcc_tree, wait=TRUE)
child.ntip <- data.frame(node.id=(1:Nnode(treeP))+n.tip,
                         n.tip=unlist(lapply(subtr.list, Ntip)))
