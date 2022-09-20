pb <- pbtree(n = 500, nsim = 1)
pbtrees <- rep(pb, 500)
for (i in 1:500){
  pbtrees[[i]] <- ape::keep.tip(pbtrees[[i]], 
                                tip = sample(pbtrees[[i]]$tip.label,
                                             100))
}
ltt95_pbtrees <- ltt95(pbtrees)
plot(ltt95_pbtrees, log = T, method = "lineages", shaded = T)
