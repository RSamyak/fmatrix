library(phylodyn)
library(ape)
library(latex2exp)
library(tidyverse)

path <- "./R/"

for(filename in list.files(path)) source(file.path(path, filename))


####Example
n <- 8
res<-construct_trees(n)
nrow(res$res)


#Indeed for n=5, there are 5 trees
#we construct the Fmatrices

F.list<-list()
for (j in 1:nrow(res$res)){
  F.list[[j]]<- Fmat_from_construct(n, res$res[j,], res$F.matrixRel)
}

frechet_var_sample(F.list)

#####################

probs <- unlist(lapply(F.list, kingman_likelihood))

probs2 <- unlist(lapply(F.list, beta_BF_likelihood, betas = 0))

frechet_var_pop(F.list, probs = probs)



betas <- c(seq(-1, 0, by = .05), 1, 2, 5, 10, 100, 1e3)

out <- rep(NA, length(betas))
out_meandist <- rep(NA, length(betas))
out_fmean_list <- list()

for(i in 1:length(betas)){
  probs <- unlist(lapply(F.list, beta_BF_likelihood, betas = betas[i]))
  
  current_out <- frechet_var_pop(F.list, probs = probs)
  
  out[i] <- current_out$var
  out_meandist[i] <- current_out$meandist
  out_fmean_list[[i]] <- current_out$fmean
}



plotF.list.png(out_fmean_list, folder = "new1")

# mylog <- function(u){
#   sign(u) * log(1 + abs(u))
# }

p1 <- ggplot(data = NULL, aes(x = as.factor(betas), y = out)) + 
  geom_point() + 
  geom_path(aes(group = 1)) +
  theme_bw() +
  scale_x_discrete(breaks = c(-1, -.75, -.5, -.25, 0, 1, 2, 5, 10, 100, 1e3)) +
  labs(x = TeX("beta"),
       y = "Frechet Variance")

p1
ggsave("figs/frechet_var_bfbeta.png", plot = p1)


p1b <- ggplot(data = NULL, aes(x = as.factor(betas), y = out_meandist)) + 
  geom_point() + 
  geom_path(aes(group = 1)) +
  theme_bw() +
  scale_x_discrete(breaks = c(-1, -.75, -.5, -.25, 0, 1, 2, 5, 10, 100, 1e3)) +
  labs(x = TeX("beta"),
       y = "Average Distance to Frechet Mean")

p1b
ggsave("figs/frechet_meandist_bfbeta.png", plot = p1b)


## does not work
# for(i in 1:length(out_fmean_list)){
#   p1 <- p1 + annotation_raster(readPNG(paste0("png/file", i, ".png")), xmin = betas[i]-.05, xmax = betas[i]+.05, ymin = out[i]-1, ymax = out[i] - 0.5)
# }



################################


load("data/F-lists.RData")


ns <- 5:9

out2 <- out2_meandist <- rep(NA, length(ns))
out2_list <- list()

for(i in 1:length(ns)){
  probs <- unlist(lapply(get(paste0("F.list",ns[i])), kingman_likelihood))
  
  current_out <- frechet_var_pop(get(paste0("F.list",ns[i])), probs = probs)
  out2[i] <- current_out$var
  
  out2_meandist[i] <- current_out$meandist
  
  out2_list[[i]] <- current_out$fmean
}


# probs <- unlist(lapply(get(paste0("F.list",10)), kingman_likelihood))
# out2_10 <- frechet_var_pop(get(paste0("F.list",10)), probs = probs)


p2 <- ggplot(data = NULL) + 
  geom_point(aes(x = as.factor(ns), y = out2)) + 
  theme_bw() +
  scale_x_discrete(breaks = ns) + 
  labs(x = TeX("n"),
       y = "Frechet Variance") + 
  ggtitle("Kingman Model")

p2

ggsave("figs/frechet_var_kingman_n.png", plot = p2)

p2b <- ggplot(data = NULL) + 
  geom_point(aes(x = as.factor(ns), y = out2_meandist)) + 
  theme_bw() +
  scale_x_discrete(breaks = ns) + 
  labs(x = TeX("n"),
       y = "Average Distance to Frechet Mean") + 
  ggtitle("Kingman Model")

p2b

ggsave("figs/frechet_meandist_kingman_n.png", plot = p2b)



out3 <- tibble()

for(i in 1:length(ns)){
  
  myFlist <- get(paste0("F.list", ns[i]))
  
  betas <- get("betas")
  
  out3_n <- tibble()
 
  for(j in 1:length(betas)){
    probs <- unlist(lapply(myFlist, beta_BF_likelihood, betas = betas[j]))
    
    current_out <- frechet_var_pop(myFlist, probs = probs)
    
    out3_n <- bind_rows(out3_n,
                        tibble(n = ns[i],
                               beta = betas[j],
                               var = current_out$var,
                               meandist = current_out$meandist,
                               fmean = list(current_out$fmean)
                               )
                        )
  }
  
  plotF.list.png(out3_n$fmean, folder = paste0("myplot2_", ns[i]))
  
  out3 <- bind_rows(out3, out3_n)
  
}


p3 <- ggplot(out3, aes(x = as.factor(round(beta, 2)), y = var, colour = as.factor(n), group = as.factor(n))) + 
  theme_bw() + 
  geom_point()+
  geom_path() + 
  scale_x_discrete(breaks = c(-1, -.75, -.5, -.25, 0, 1, 2, 5, 10, 100, 1e3))  +
  # scale_x_log10() + 
  labs(x = TeX("beta"),
       y = "Frechet Variance",
       colour = "n")

p3
ggsave("figs/frechet_var_bfbeta_n.png", plot = p3)


p3b <- ggplot(out3, aes(x = as.factor(round(beta, 2)), y = meandist, colour = as.factor(n), group = as.factor(n))) + 
  theme_bw() + 
  geom_point()+
  geom_path() + 
  scale_x_discrete(breaks = c(-1, -.75, -.5, -.25, 0, 1, 2, 5, 10, 100, 1e3))  +
  # scale_x_log10() + 
  labs(x = TeX("beta"),
       y = "Average Distance to Frechet Mean",
       colour = "n")

p3b
ggsave("figs/frechet_meandist_bfbeta_n.png", plot = p3b)



# p4 <- ggplot(out3) + 
#   theme_bw() + 
#   geom_point(aes(x = as.factor(n), y = var, colour = as.factor(beta)))+
#   scale_x_discrete(breaks = ns)  + 
#   labs(x = TeX("n"),
#        y = "Frechet Variance",
#        colour = "beta")
# 
# p4
# ggsave("figs/frechet_var_n_bfbeta.png", plot = p4)
