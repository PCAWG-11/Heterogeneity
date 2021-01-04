#### looking at TRACERx 100 trees and extracting averages

library(data.tree)
library(ggplot2)

# read in Tx100 trees and mutation burdens
treerda <- "Tx100TreesCRUK.rda"
load(file = treerda)

totalmutburden <- read.table(file = "Tx100MutbBurden.tsv", as.is = T)
totalmutburden <- setNames(object = totalmutburden$mutburden, nm = rownames(totalmutburden))


# simulation
nsamples <- length(listTrees)
nreps <- 1000
nsims <- nsamples*nreps


ischildv <- rep(NA, nsims)
for (i in 1:nsims) {
  
  # take a tree
  singletree <- listTrees[[i %% nsamples + 1]]
  # print(singletree$levelName)
  
  # get the number of subclones, skip if < 2
  subclones <- singletree$Get(attribute = "name")[-c(1,2)]
  if (length(subclones) < 2) next
  
  # randomly pick 2 of the subclones 
  subsub <- as.integer(sample(x = subclones, size = 2, replace = F))
  
  # and check whether one is the child of the other (i.e. linear)
  ischildv[i] <- !is.null(FindNode(FindNode(node = singletree, name = subsub[1]), name = subsub[2])) |
    !is.null(FindNode(FindNode(node = singletree, name = subsub[2]), name = subsub[1]))
  
}


# get the odds ratio branching/linear
out <- c(by(data = ischildv, INDICES = rep(1:nreps, rep(nsamples,nreps)), FUN = function(x) sum(!x, na.rm = T)/sum(x, na.rm = T)))
print(quantile(x = out, probs = c(0.025, .5, 0.975)))

# get the frequencies branching, linear
out2 <- c(by(data = ischildv, INDICES = rep(1:nreps, rep(nsamples,nreps)), FUN = function(x) sum(!x, na.rm = T)/length(x)))
txprobs <- quantile(x = out2, probs = c(0.025, .5, 0.975))

# plot 
plottx <- data.frame(type = c("branching", "linear"), lo = c(txprobs[[1]], 1-txprobs[[3]]), med = c(txprobs[[2]], 1-txprobs[[2]]), hi = c(txprobs[[3]],1-txprobs[[1]]))
pa <- ggplot(data = plottx, mapping = aes(x = type)) + geom_col(mapping = aes(y = med)) + geom_errorbar(mapping = aes(ymin = lo, ymax = hi), width = .25)
pa <- pa + theme_minimal() + theme(axis.title.x = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.y = element_blank()) + labs(y = "Probability in TRACERx 100") + ylim(c(0,38/37))
# pa
ggsave(filename = "figure_4e.pdf", plot = pa, width = 2, height = 4)


### check relation linear/branching vs mutation burden
lindf <- data.frame(prob_linear = c(by(data = ischildv, INDICES = rep(names(listTrees)[c(2:nsamples,1)], nreps), sum)/nreps))
lindf$crukid <- row.names(lindf)
lindf$mutburden <- c(totalmutburden[lindf$crukid])

corout <- cor.test(x = lindf$prob_linear, y = lindf$mutburden, use = "na.or.complete", method = "spearman")

p1 <- ggplot(data = lindf[!is.na(lindf$prob_linear), ], mapping = aes(x = mutburden, y = prob_linear)) + geom_vline(xintercept = 10^c(2:3), size = .2, color = "gray") +
  geom_point(stroke = 0) + annotation_logticks(sides = "b") + lims(y = c(0,1)) + scale_x_log10(lim = c(50,5000), breaks = 10^c(2:5)) + 
  annotate(geom = "text", x = 4000, y = 0.95, label = paste0("\u03C1 = ", round(corout$estimate, digits = 3), "\n", "P = ", round(corout$p.value, digits = 3))) +
  theme_minimal() + theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank()) + labs(x = "Mutation burden (exome)", y = "Proportion linear relationships")
# p1
ggsave(filename = "figure_4g.pdf", plot = p1, width = 6, height = 4)
