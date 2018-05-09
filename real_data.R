effects <- data.frame(Group = factor(c(rep("Age", 2), 
                                          rep("Sex", 2), 
                                          rep("Diabetes", 2), 
                                          rep("Stroke or TIA", 2),
                                          rep("Creatinine", 3), 
                                          rep("CHADS2", 3), 
                                          rep("VKA", 2), 
                                          rep("TTR", 2)),
                                     levels = c("Age", "Sex", "Diabetes", "Stroke or TIA",
                                                "Creatinine", "CHADS2",
                                                "VKA", "TTR"),
                                        ordered = T),
                      Subgroup = c(1:2, 1:2, 1:2, 1:2,
                                   1:3, 1:3, 1:2, 1:2),
                      pos.NOAC = c(496, 415, 382, 531,
                                   622, 287, 483, 428,
                                   249, 405, 256, 
                                   69, 247, 596,
                                   386, 522, 509, 313),
                      total.NOAC = c(18073, 11188, 10941, 18371,
                                     20216, 9096, 20699, 8663,
                                     5539, 13055, 10626,
                                     5058, 9563, 14690,
                                     13789, 15514, 16219, 12742),
                      pos.war = c(578, 532, 478, 634,
                                  755, 356, 615, 495,
                                  311, 546, 255,
                                  90, 290, 733,
                                  513, 597, 653, 392),
                      total.war = c(18004, 11095, 10839, 18390,
                                    20238, 8990, 20637, 8635,
                                    5503, 13155, 10533,
                                    4942, 9757, 14528,
                                    13834, 15395, 16297, 12904))

ORpvalues <- function(data) {
  p.vec <- apply(data[, -(1:2)], 1, function(v){
                        v <- as.numeric(v)
                        v[c(2, 4)] <- v[c(2, 4)] - v[c(1, 3)]
                        temp <- matrix(v, nrow = 2)
                        temp1 <- fisher.test(temp)
                        return(c(temp1$estimate,
                                 temp1$p.value))
                      })
  return(p.vec)
}

## compute the ods ratio and p-value in table Table 2(a)
p.vec <- ORpvalues(effects)



effects <- data.frame(effects, odds.ratio = p.vec[1, ], pvalue = p.vec[2, ])

effect.tex <- effects

effect.tex$Subgroup <- c("leq 75", "geq 75",
                      "Female", "Male", 
                      "No", "Yes",
                      "No", "Yes",
                      "leq 50", "50-80", "geq 80",
                      "0-1", "2", "3-6",
                      "Naive", "Experienced",
                      "leq 66", "geq 66")
effect.tex$NOAC <- apply(effects[, 3:4], 1, paste, collapse = "/")
effect.tex$war <- apply(effects[, 5:6], 1, paste, collapse = "/")
effect.tex <- effect.tex[, c(2, 9:10, 7:8)] 


## compute Fisher BHPC p-value for r = 1, 2, 3 within each group ##########
## p.grp generates Table 2(c) ####

gr.idx <- levels(effects$Group)

p.grp <- sapply(1:length(gr.idx), function(i){
                        v <- effects$pvalue[effects$Group == gr.idx[i]]
                        v.p1 <- pchisq(-2 * sum(log(v)), 2 * length(v), lower.tail = F)
                        v1 <- v[-which.min(v)]
                        v.p2 <- pchisq(-2 * sum(log(v1)), 2 * length(v1),
                                       lower.tail = F)
                        return(c(v.p1, v.p2))
                      })
p.grp <- rbind(p.grp, rep(1, 8))
p.grp[3, 5] <- max(effects$pvalue[effects$Group == gr.idx[5]])
p.grp[3, 6] <- max(effects$pvalue[effects$Group == gr.idx[6]])
colnames(p.grp) <- gr.idx


##### calculate the new GBHPC p-values #######
p.grp <- p.grp[, sort(p.grp[1, ], index.return = T)$ix]
p.grp1 <- p.grp[, sort(p.grp[2, ], index.return = T)$ix]



p.Bon <- sort(effects$pvalue)[-1] * (17:1)
l2 <- sort(p.grp[2, ]) * 8:1
l2 <- sapply(1:length(l2), function(i) max(l2[max((i-2), 1):i]))
p.new <- c(p.grp[1, -1]* 8,
           l2,
           sort(p.grp1[3, ])[1:2] * c(2, 1))
p.mat <- cbind(p.Bon, p.new)
rownames(p.mat) <- 2:18
require(xtable)
print(xtable(p.mat, digits = 2, display = c("s", "E", "E")))



