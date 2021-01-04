### Code for Dentro et al. Figure 4A-D and F

library(ggplot2)
library(reshape2)
library(boot)
# library(gridExtra)

source("pcawg.colour.palette.R")


linbranchresults_rep <- read.delim(file = "additionalITH_phasingresults.tsv", as.is = T)

# PCAWG colour scheme
cvect <- pcawg.colour.palette(x = tolower(sub(pattern = "-", replacement = ".", x = unique(linbranchresults_rep$histology_abbreviation))), scheme = "tumour.subtype")
names(cvect) <- unique(linbranchresults_rep$histology_abbreviation)
cvect[c("Skin-Melanoma-Acral", "Skin-Melanoma-Cut", "Skin-Melanoma-Mucosal", "Kidney-RCC-Clear", "Kidney-RCC-Pap")] <- c("#000000", "#000000", "#000000", '#FF4500', '#FF4500')


# manuscript checks and numbers
print("tumors with sufficient power and at least one phaseable pair in the correct context for linear phasing")
sum(linbranchresults_rep$callable_lin > 0)

print("show discordant in-cis SNVs pairs")
sum(linbranchresults_rep$lin_total > 0)
print("clone-subclone")
sum(linbranchresults_rep$lin_cl_subcl > 0)

print("subclone-other subclone")
sum(linbranchresults_rep$lin_diff_subcl > 0)
print("subclone-same subclone")
sum(linbranchresults_rep$lin_same_subcl > 0)


print("tumors with sufficient power and at least one phaseable pair in the correct context")
sum(linbranchresults_rep$callable_branching > 0)

print("show discordant in-trans SNVs pairs")
sum(linbranchresults_rep$bra_total > 0)

print("subclone-other subclone")
sum(linbranchresults_rep$bra_diff > 0)
print("subclone-same subclone")
sum(linbranchresults_rep$bra_same > 0)

print("subclone-subclone both")
sum(linbranchresults_rep$bra_diff > 0 & linbranchresults_rep$bra_same > 0)




# test to see whether the number of subclones called is associated with the total number of linear & branching pairs
outtab <- data.frame(subrecon = linbranchresults_rep$num_subclones,
                     phasing = linbranchresults_rep$lin_total + linbranchresults_rep$bra_total > 0, stringsAsFactors = F)

table(outtab$subrecon, outtab$phasing)
chisqout <- chisq.test(table(outtab$subrecon, outtab$phasing)[-5,])

p1 <- ggplot() + geom_bar(data = outtab, mapping = aes(fill = phasing, x = subrecon), position = "fill")
p1 <- p1 + geom_text(data = data.frame(table(outtab$subrecon)), mapping = aes(x = as.numeric(Var1)-1, y = 1, label = Freq), nudge_y = .05)
p1 <- p1 + scale_y_continuous(labels = scales::percent, breaks = (0:4)/4) + labs(x = "# Consensus subclones", y = "Fraction of samples") +
  scale_fill_brewer(type = "qual") + theme_minimal() +
  guides(fill = guide_legend("Has ≥ 1 linear or\n branching pair"))
p1 <- p1 + annotate(geom = "text", x = 2, y = 1.1, label = paste0("p-value = ", sprintf(fmt = "%e", chisqout$p.value)))
# p1

ggsave(filename = "figure_4c.pdf", plot = p1, width = 6, height = 6)





## visualise fraction of callable samples which has evidence for each type of linear evo
callabletypes <- sort(unique(linbranchresults_rep[linbranchresults_rep$lin_total >= 1 | linbranchresults_rep$bra_total >= 1, "histology_abbreviation"]))

plotdf1 <- linbranchresults_rep[linbranchresults_rep$lin_total >= 1, ]
plotdf1 <- melt(data = plotdf1, id.vars = c("samplename", "histology_abbreviation"), measure.vars = c("lin_cl_subcl", "lin_same_subcl", "lin_diff_subcl"))
plotdf1 <- plotdf1[plotdf1$value > 0, ]

# downweight samples if they've got multiple types of linear pairs so that sum of the slices == 1
plotdf1$weight <- c(1/table(plotdf1$samplename)[plotdf1$samplename])


# downweight samples so that total height (sum of all callable samples) == 1
plotdf2 <- data.frame(table(linbranchresults_rep[linbranchresults_rep$callable_lin >= 1 &
                                                   linbranchresults_rep$histology_abbreviation %in% callabletypes, "histology_abbreviation"]))
rownames(plotdf2) <- plotdf2$Var1
plotdf1$weight <- plotdf1$weight/plotdf2[plotdf1$histology_abbreviation, "Freq"]


plotdf1$variable <- factor(x = plotdf1$variable, levels = c("lin_diff_subcl", "lin_same_subcl", "lin_cl_subcl"))

# consistent order
plotdf1$histology_abbreviation <- factor(plotdf1$histology_abbreviation, levels = callabletypes)
plotdf2$histology_abbreviation <- factor(plotdf2$Var1, levels = callabletypes)

p1 <- ggplot(data = plotdf1, mapping = aes(x = histology_abbreviation)) + geom_bar(mapping = aes(weight = weight, fill = variable, alpha = variable), show.legend = F)
p1 <- p1 + geom_text(data = plotdf2, mapping = aes(x = Var1, y = 1, label = Freq), nudge_y = .05)
p1 <- p1 + scale_x_discrete(drop = F)
p1 <- p1 + scale_fill_manual(values = c(callable = '#ffffff', lin_diff_subcl = '#fbb4ae', lin_same_subcl = '#b3cde3', lin_cl_subcl = '#bee7b5'))
p1 <- p1 + scale_alpha_manual(values = c(callable = 0, lin_diff_subcl = 1, lin_same_subcl = 1, lin_cl_subcl = 1))
# p1 <- p1 + scale_fill_manual(values = c(callable = '#ffffff', lin_diff_subcl = '#fbb4ae', lin_same_subcl = '#b3cde3', lin_cl_cl = '#ccebc5', lin_cl_cl_xor = '#decbe4'))
p1 <- p1 + theme_minimal() + theme(axis.text.x = element_text(angle = 90), axis.title.x = element_blank(), panel.grid.minor.y = element_blank()) + labs(y = "fraction of samples")
# p1

ggsave(filename = "figure_4a.pdf", plot = p1, width = 9, height = 3)




## visualise fraction of callable samples which has evidence for each type of branching evo
plotdf1 <- linbranchresults_rep[linbranchresults_rep$bra_total >= 1, ]
plotdf1 <- melt(data = plotdf1, id.vars = c("samplename", "histology_abbreviation"), measure.vars = c("bra_same", "bra_diff"))
plotdf1 <- plotdf1[plotdf1$value > 0, ]
plotdf1$weight <- c(1/table(plotdf1$samplename)[plotdf1$samplename])

plotdf2 <- data.frame(table(linbranchresults_rep[linbranchresults_rep$callable_branching >= 1 &
                                                   linbranchresults_rep$histology_abbreviation %in% callabletypes, "histology_abbreviation"]))
rownames(plotdf2) <- plotdf2$Var1
plotdf1$weight <- plotdf1$weight/plotdf2[plotdf1$histology_abbreviation, "Freq"]

plotdf1$variable <- factor(x = plotdf1$variable, levels = c("bra_same", "bra_diff"))

plotdf1$histology_abbreviation <- factor(plotdf1$histology_abbreviation, levels = callabletypes)
plotdf2$histology_abbreviation <- factor(plotdf2$Var1, levels = callabletypes)

p1 <- ggplot(data = plotdf1, mapping = aes(x = histology_abbreviation)) + geom_bar(mapping = aes(weight = weight, fill = variable, alpha = variable), show.legend = F)
p1 <- p1 + geom_text(data = plotdf2, mapping = aes(x = Var1, y = 1, label = Freq), nudge_y = .05)
p1 <- p1 + scale_x_discrete(drop = F)
# p1 <- p1 + scale_fill_manual(values = c(callable = '#ffffff', bra_same_cl_cl = '#fbb4ae', bra_same_sub_sub = '#b3cde3', bra_diff_cl_sub = '#ccebc5', bra_diff_sub_sub = '#decbe4'))
p1 <- p1 + scale_fill_manual(values = c(callable = '#ffffff', bra_diff = '#fbb4ae', bra_same = '#b3cde3'))
p1 <- p1 + scale_alpha_manual(values = c(callable = 0, bra_diff = 1, bra_same = 1))
p1 <- p1 + theme_minimal() + theme(axis.text.x = element_text(angle = 90), axis.title.x = element_blank(), panel.grid.minor.y = element_blank()) + labs(y = "fraction of samples")
# p1

ggsave(filename = "figure_4b.pdf", plot = p1, width = 9, height = 3, useDingbats=FALSE)




### check lin vs branching probabilities: even the field by disallowing assignments to subclones and check in callable branching areas only
plotdf1 <- linbranchresults_rep[linbranchresults_rep$callable_sub_sub_bra > 0, ]
plotdf1$call <- ifelse(plotdf1$lin_subsub_bra > 0, ifelse(plotdf1$bra_subsub_bra > 0, "both", "linear"), ifelse(plotdf1$bra_subsub_bra > 0, "branching", "callable"))

table(plotdf1$call)
wilcout <- wilcox.test(x = rowSums(plotdf1[plotdf1$call == "branching", c("num_clonal", "num_subclonal")]), y = rowSums(plotdf1[plotdf1$call == "linear", c("num_clonal", "num_subclonal")]))

plotdf1$call <- factor(x = plotdf1$call, levels = c("callable", "branching", "linear", "both"))



# check effect of mutation burden on lin vs branching
p1 <- ggplot(data = plotdf1, mapping = aes(x = call, y = num_clonal + num_subclonal)) + geom_boxplot(outlier.shape = NA, width = .25, fill = "gray90") + geom_jitter(stroke = 0, width = .05, alpha = .35) + scale_y_log10() + annotation_logticks(sides = "l")
p1 <- p1 + theme_minimal() + theme(axis.text.x = element_text(angle = 45), panel.grid.minor.y = element_blank()) +  labs(x = "", y = "Mutation burden")
p1 <- p1 + annotate(geom = "segment", x = c(2,2,3), y = c(4e6,3.5e6, 3.5e6), xend = c(3,2,3), yend = c(4e6,4e6, 4e6)) + annotate(geom = "text", x = 2.5, y = 5.5e6, label = paste0("P = ", round(wilcout$p.value, 3)))
# p1
ggsave(filename = "figure_4f.pdf", plot = p1, width = 5, height = 4)



nbra_boot <- boot(data = plotdf1, statistic = function(x, i) {c(sum(x[i,"call"] %in% c("linear", "both")),
                                                                sum(x[i,"call"] %in% c("branching", "both")))}, R = 100000, parallel = "multicore", ncpus = 4)

bootdf <- data.frame(call = c("linear", "branching"), t0 = nbra_boot$t0)
bootdf[, c("ymin", "ymax")] <- t(apply(X = nbra_boot$t, MARGIN = 2, FUN = quantile, probs = c(.025, .975)))
quantile(nbra_boot$t[,1]/nbra_boot$t[,2], probs = c(.025, .975))
nbra_boot$t0[1]/nbra_boot$t0[2]


# duplicate rows with both
bothlinear <- plotdf1[plotdf1$call == "both",]
bothbranching <- bothlinear
bothlinear$call <- "linear"
bothbranching$call <- "branching"
plotdf1 <- rbind(plotdf1[plotdf1$call != "both", ], bothlinear, bothbranching)

p1 <- ggplot(data = plotdf1[plotdf1$call != "callable", ]) + geom_bar(mapping = aes(x = call, fill = histology_abbreviation, color = samplename %in% bothlinear$samplename, size = samplename %in% bothlinear$samplename), show.legend = F)
p1 <- p1 + geom_errorbar(data = bootdf, mapping = aes(x = call, ymin = ymin, ymax = ymax), width = .25) + scale_size_manual(values = c('TRUE' = 1, 'FALSE' = 0)) + scale_color_manual(values = c('TRUE' = "orange", 'FALSE' = "white"))
p1 <- p1 + scale_fill_manual(values = cvect) + theme_minimal() + theme(axis.title.x = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.y = element_blank()) + labs(y = "# Samples") + ylim(c(0,38))
# p1

ggsave(filename = "figure_4d.pdf", plot = p1, width = 2, height = 4)


