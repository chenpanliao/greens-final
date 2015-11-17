setwd("/Users/apan/Dropbox/greens-final")
require(PMCMR)
require(caret)
require(relaimpo)
require(QuantPsyc)
require(vegan)
require(xtable)
require(magrittr)
require(tidyr)
require(ggplot2)
require(hydroGOF)
require(car)
library(MuMIn)
mds.try.n <- 20
NRMSE <- function(sim, obs){
  hydroGOF::rmse(sim, obs) / sd(obs)
}
multimerge <- function(x, ...) {
  if(!is.list(x)) stop("x must be a list.")
  y <- x[[1]]
  i <- 2
  while(i <= length(x)){
    y <- merge(y, x[[i]], by='row.names', all=T, sort=T)
    rownames(y) <- y$Row.names
    y$`Row.names` <- NULL
    colnames(y) <- names(x)[1:i]
    i <- i + 1
  }
  return(y)
}
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  ra <- abs(cor(x, y, method="kendall"))
  r <- cor(x, y, method="kendall")
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 1.5/strwidth(txt)
  text(0.5, 0.5, txt, cex = ra*1.2+0.5)
}
panel.smooth <- function (x, y, col = par("col"), bg = NA, pch = par("pch"),
  cex = 1, col.smooth = "red", span = 4/5, iter = 100, ...)
{
  points(x, y, pch = pch, col = col, bg = bg, cex = cex)
  ok <- is.finite(x) & is.finite(y)
  if (any(ok))
    lines(stats::lowess(x[ok], y[ok], f = span, iter = iter),
      col = col.smooth, ...)
}

quartzFonts(sans = quartzFont(rep("Noto Sans CJK TC Regular", 4)),
            serif = quartzFont(rep("Noto Sans CJK TC Regular", 4)))
BAF.mat <- matrix(
  c(1.4, 4.4, 3.9, 4.9, 4.3, 4.5, 3.7, 6.1, 7, 6.4, 4.3, 3.4, 5.4), 1,
  dimnames = list(c("BAF"), letters[1:13]) )

d <- read.csv("dat/All_Data_20151115.csv")
names(d) %<>% gsub("^G([a-z])$", "\\1", .)
names(d)[grep("^[a-z]$", names(d))] %<>% toupper
names(d)[grep("^[a-z]$", names(d))] %<>% toupper
# bio data 為 NA 者全填為 0，但以後要再檢查一次
d <- d[, c(1:17, 19,21,23,25, 18,20,22,24)]
d[, grep("^(Dgs|Hs)(Spider|Insect|Grass|Wood)$", names(d))] [
  is.na(d[, grep("^(Dgs|Hs)(Spider|Insect|Grass|Wood)$", names(d))])
] <- 0


# 所有樣本
# 有 NA
dG0 <- d[, grep("^[A-X]$", names(d))]
dG1 <- d[, grep("^[A-W]$", names(d))]
dB0 <- d[, grep("^(Dgs|Hs)[a-zA-Z]+$", names(d))]


## 找出總面積不正常的點
quartz(width=7, height=3)
par(cex=10/12, mar=c(4,4,0,0)+0.1, family = "Noto Sans CJK TC")
hist(rowSums(dG0), nclass=40, main="", xlim=c(0,60000), ylim=c(0,500), col=3, xlab="總面積（m²）", ylab="次數")
  abline(v=c(8000,12000), lty=2, col=2)
  text(4000, 260, "< 8000\n(219/1169)")
  text(16000, 260, "> 12000\n(26/1169)")
# quartz.save("./slide/invalid-area.pdf", type="pdf")
quartz.save("./slide/invalid-area.png", dpi=300, type="png")
dev.off()


## 找出未定義面積不正常的點
quartz(width=7, height=3)
par(cex=10/12, mar=c(4,4,0,0)+0.1, family = "Noto Sans CJK TC")
hist(dG0$X, nclass=40, main="", xlim=c(0,20000), ylim=c(0,600), col=3, xlab="未定義面積（m²）", ylab="次數")
  abline(v=c(8000), lty=2, col=2)
  text(12000, 550, "> 8000\n(125/1169)")
# quartz.save("./slide/invalid-gx-area.pdf", type="pdf")
quartz.save("./slide/invalid-gx-area.png", dpi=300, type="png")
dev.off()


## 有問題的景觀資料
Location <- as.character(d$Location)
  Location[duplicated(Location)] <- ""
quartz(width=7, height=7)
par(family = "Noto Sans CJK TC", cex=10/12)
hv <- heatmap(
  log10(dG0 + 1) %>% as.matrix,
  Rowv = NA, Colv = NA,
  # col = rainbow(256),
  col = gray(seq(1, 0, length=6)),
  scale = "none", na.rm = F,
  margins = c(6,6),cexRow=0.3, cexCol=1,
  labRow = Location,
  xlab = "", ylab =  "")
# quartz.save("./slide/invalid-type.pdf", type="pdf")
quartz.save("./slide/invalid-type.png", dpi=300, type="png")
dev.off()


## 有問題的多樣性資料
Location <- as.character(d$Location)
  Location[duplicated(Location)] <- ""
quartz(width=7, height=7)
par(family = "Noto Sans CJK TC", cex=10/12)
hv <- heatmap(
  log10(dB0 + 1) %>% as.matrix ,
  Rowv = NA, Colv = NA,
  # col = terrain.colors(50),
  col = gray(seq(1, 0, length=6)),
  scale = "none", na.rm = F,
  margins = c(6,6),cexRow=0.3,cexCol=1,
  labRow = Location,
  labCol = paste(c("木本","草本","蜘蛛","昆蟲"), c(rep("Simpson",4), rep("Shannon",4)), sep="\n"),
  xlab = "", ylab =  "")
# quartz.save("./slide/invalid-bio.pdf", type="pdf")
quartz.save("./slide/invalid-bio.png", dpi=300, type="png")
dev.off()


## Simpson 超過 1 不合理
Location <- as.character(d$Location)
  Location[duplicated(Location)] <- ""
quartz(width=7, height=7)
par(family = "Noto Sans CJK TC", cex=10/12)
hv <- heatmap(
  ((dB0[, grep("Dgs[a-zA-z]+", names(dB0))] > 1)*2) %>% as.matrix ,
  Rowv = NA, Colv = NA,
  # col = terrain.colors(50),
  col = gray(seq(1, 0, length=2)),
  scale = "none", na.rm = F,
  margins = c(6,6),cexRow=0.3,cexCol=1,
  labRow = Location,
  labCol = paste(c("木本","草本","蜘蛛","昆蟲"), rep("Simpson",4), sep="\n"),
  xlab = "", ylab =  "")
# quartz.save("./slide/invalid-bio.pdf", type="pdf")
quartz.save("./slide/invalid-bio-simpson.png", dpi=300, type="png")
dev.off()












# 求解用的 dataset
# 無 NA
# 多樣性資料只看 shannon
# 去除 Gl
# 朝陽4 NO101 是唯一 K != 0 的點
# 去除 < 8000 | > 12000 area
# 去除 Gx > 8000
# 面積以列標準化為 sum = 1
# NO=="640" 很怪
# isAreaNormalization <- F
isAreaMeanCentering <- T
isTotalSameArea <- T # 總和調成 1
isUnknownAreaIgnored <- T # 不理 X 其它總和調成 1
dt <- d %>%
  na.omit %>%
  subset(., X+A+B+C+D+E+F+G+H+I+J+K+L+M >= 5000) %>%
  subset(., X+A+B+C+D+E+F+G+H+I+J+K+L+M <= 15000) %>%
  subset(., X <= 5000)
dim(dt); str(dt)
# 環境資料
# 一次項
dG1t1 <- dt %>% .[, grep("^[A-Z]$", names(.)) ]
if (isUnknownAreaIgnored == T)
  dG1t1 %<>% .[, !colnames(.) %in% "X" ]
if (isTotalSameArea)
  dG1t1 %<>% decostand(., "total", 1)
dG1t1Mean <- dG1t1 %>% colMeans # 求 col mean
if (isAreaMeanCentering)
  dG1t1 %<>% apply(., 2, function(x){x-mean(x)}) # 中心化
dG1t1 %<>% .[, !colnames(.) %in% c("X", "G", "H", "K", "L", "M")  ] %>%
  as.data.frame
  dim(dG1t1); head(dG1t1)
# # 二次項
dG1t2 <- dG1t1^2 %>% as.data.frame
  names(dG1t2) %<>% paste0(., "2")
  dim(dG1t2); head(dG1t2)
# 交互作用項
dG1tI <- model.matrix(~(.)^2, dG1t1) %>%
  as.data.frame %>%
  .[, grep(":", colnames(.))] %>%
  as.data.frame
names(dG1tI) %<>% gsub(":", "", .) ; dim(dG1tI); head(dG1tI)
# 結合三者 # 減去沒用的交互作用
dG1t <- cbind(dG1t1, dG1t2, dG1tI) %>%
  as.data.frame %>%
  .[, colSums(.) != 0] %>% .[!colnames(.) %in% c("A2","B2","C2","D2","E2","F2","I2","J2","AC","AF","BC","BF","CD","CE","CF","CI","CJ","DF","EF","EI","EJ","FI","FJ")]
str(dG1t)

# 生態資料
# 只用 Shaanon
# 木本拿掉
dB0t <- dt %>%
  .[, grep("^Hs(Spider|Insect|Grass)$", names(.))] %>%
  decostand(., "range", 2)
str(dB0t)

# 所有採用資料圖
quartz(width=7, height=3.5)
par(cex=8/12, mar=c(5,5,0,0)+0.1, family = "Noto Sans CJK TC", mgp=c(4,1,0))
boxplot(cbind(dB0t, dG1t), xlab="變數", ylab="標準化後數值", las=2)
  # axis(1, bp, names(dG1t))
# quartz.save("./slide/var-summary.pdf", type="pdf")
quartz.save("./slide/var-summary.png", dpi=300, type="png")
dev.off()





## 求原始生態資料降維後分數 Bscore.o
stepDirection <- "both"
# stepDirection <- "forward"
Bscore.o <- data.frame(NO = dt$NO)
fitPCA <- function(dt){
  fit <- princomp(dt, cor=F)
  score <- fit$scores[,1]
  loading <- t(loadings(fit)[,1])
  if(max(loading) < 0) {score <- -score; loading <- -loading}
  print(summary(fit))
  print(round(loading, digits=6))
  return((score - mean(score))/sd(score))
}
fitFA <- function(dt){
  fit <- factanal(dt, 1, rotation="varimax", scores="regression", trace=T)
  score <- as.vector(fit$scores)
  loading <- t(loadings(fit))
  if(max(loading) < 0) {score <- -score; loading <- -loading}
  fit
  print(t(loadings(fit)), digits=6)
  return((score - mean(score))/sd(score))
}
fitMDS <- function(dt){
  fit <- vegan::metaMDS(dt, distance="euclidean", k = 1, trymax = mds.try.n)
  score <- scores(fit)[,1]
  loading <- t(fit$species[,1])
  if(max(loading) < 0) {score <- -score; loading <- -loading}
  print(round(fit$stress, digits=6))
  print(round(loading, digits=6))
  return((score - mean(score))/sd(score))
}
fitRDA <- function(dtB, dtG){
  # fit.full <- vegan::rda(dtB ~ . , data=dtG, scale=F)
  # fit.null <- vegan::rda(dtB ~ 1 , data=dtG, scale=F)
  # fit.start <- vegan::rda(dtB ~ A+B+C+D+E+F+I+J, data=dtG, scale=F)
  # fit <- ordiR2step(fit.start, scope = list(upper=formula(fit.full), lower=formula(fit.null)), direction=stepDirection)
  fit <- vegan::rda(dtB ~ ., data=dtG, scale=F)
  score <- summary(fit)$sites[,1]
  loading <- t(summary(fit)$species[,1])
  if(max(loading) < 0) {score <- -score; loading <- -loading}
  print(head(summary(fit)))
  print(round(loading, digits=6))
  return((score - mean(score))/sd(score))
}
Bscore.o$PCA <- fitPCA(dB0t)
Bscore.o$FA <- fitFA(dB0t)
Bscore.o$MDS <- fitMDS(dB0t)
Bscore.o$RDA <- fitRDA(dB0t, dG1t)

quartz(width=7, height=7)
par(mfrow=c(1,4), cex=4/12)
pairs(
  cbind(Bscore.o[,-1], dB0t),
  lower.panel=panel.smooth,
  upper.panel=panel.cor, pch=".")
quartz.save("./slide/原始資料降維分數.png", dpi=300, type="png")
# quartz.save("./slide/原始資料降維分數.pdf", type="pdf")
dev.off()







## multiple regression: diversity ~ area^1 + area^2 + interaction(among area^1)
calc.relimp(lm(Bscore.o[[1]] ~ ., dG1t1), type=c("lmg", "car", "first"), rela=F)
fit.o <-  vector(mode = "list", 8)
names(fit.o) <- paste(
  c("PCA", "PCA", "FA", "FA", "MDS", "MDS", "RDA", "RDA"),
  rep(c("1", "F"), 4),
  sep="." )
stepDirection <- "both"
# this.k <- rep(c( log(nrow(dG1t)) , log(nrow(dG1t)) ), 4)
this.k <- rep(c(2, 2), 4)
for(i in 1:length(fit.o)) {
  p <- c(1,1,2,2,3,3,4,4)
  m1 <- lm(Bscore.o[,-1][[ p[i] ]] ~ A+B+C+D+E+F+I+J , data=dG1t)
  m0 <- lm(Bscore.o[,-1][[ p[i] ]] ~ A+D+E+F+I , data=dG1t)
  m2 <- lm(Bscore.o[,-1][[ p[i] ]] ~ . , data=dG1t)
  if(i %% 2 == 1) {
    # 從 reduced model m1 開始, AIC
    fit.o[[i]] <- stats::step(m1, scope = list(upper = formula(m2), lower = formula(m0)), direction = stepDirection, trace=1, k=this.k[i])
  } else {
    # 從 full model m2 開始, AIC
    fit.o[[i]] <- stats::step(m2, scope = list(upper = formula(m2), lower = formula(m0)), direction = stepDirection, trace=1, k=this.k[i])
  }
}

# 結合係數成data.frame: coefM.o
{coefM.o <- multimerge(lapply(fit.o, function(x){data.frame(coef(x))}))} %T>% print %>%
  xtable(., digits=3, align = c("l", rep("r", ncol(.))) ) %>%
  print(., NA.string="", booktabs=F, size="normalsize", math.style.negative=T,
    file="slide/迴歸係數.tex",
    sanitize.colnames.function = function(x){gsub("(PCA|FA|MDS|RDA).(AIC|BIC)", "\\1\\\\textsubscript{\\2}", x)  },
    sanitize.rownames.function = function(x){
        x[x=="(Intercept)"] <- "Int."
        x <- gsub("^([A-Z])$", "$ \\1 $", x)
        x <- gsub("^([A-Z])2$", "$ \\1^2 $", x)
        x <- gsub("^([A-Z])([A-Z])$", " $\\1 \\\\times \\2 $", x)
        return(x)
    },
    tabular.environment = "mytabular"
  )

# show vif table
{coefVifM.o <- fit.o %>% lapply(., vif) %>% multimerge} %T>% print %>%
  xtable(., digits=3, align = c("l", rep("r", ncol(.))) ) %>%
  print(., NA.string="", booktabs=F, size="normalsize", math.style.negative=T,
    file="slide/VIF.tex",
    sanitize.colnames.function = function(x){gsub("(PCA|FA|MDS|RDA).(AIC|BIC)", "\\1\\\\textsubscript{\\2}", x)  },
    sanitize.rownames.function = function(x){
        x[x=="(Intercept)"] <- "Int."
        x <- gsub("^([A-Z])$", "$ \\1 $", x)
        x <- gsub("^([A-Z])2$", "$ \\1^2 $", x)
        x <- gsub("^([A-Z])([A-Z])$", " $\\1 \\\\times \\2 $", x)
        return(x)
    },
    tabular.environment = "mytabular"
  )



## beta : relative importance
ri <- list()
for(i in 1:length(fit.o)) {
  ri[[i]] <- calc.relimp(fit.o[[i]], type=c("lmg"), rela=F)
}
quartz(width=7, height=4)
par(mar=c(4,3,1,0)+0.1, family = "Noto Sans CJK TC", mgp=c(3,1,0), cex=8/12, mfrow=c(2,4))
for(i in 1:length(ri)) {
  barplot( ri[[i]]@lmg, xlab="", ylab="R² 貢獻量", beside=T, ylim=c(0,0.2), las=2)
  title(names(fit.o)[i])
}
# quartz.save( "slide/beta相對貢獻量.pdf", type="pdf")
quartz.save( "slide/beta相對貢獻量.png", type="png", dpi=300)
dev.off()

# # beta : relative importance (with bootstrap 95% CI)
# riBoot <- vector("list", length(fit.o))
# names(riBoot) <- names(fit.o)
# for(i in 1:length(riBoot)) {
#   set.seed(52004800)
#   riBoot[[i]] <- boot.relimp(fit.o[[i]], b = 1000, type = c("lmg"), rela = F, rank=F, diff=F)
# }
# quartz(width=10, height=5)
# par(mar=c(4,4,1,0)+0.1, family = "Noto Sans CJK TC", mgp=c(3,1,0), cex=10/12, mfrow=c(2,4))
# for(i in 1:length(riBoot)) {
#   relimpBootResult <- rbind(
#     booteval.relimp(riBoot[[i]])@lmg.upper,
#     booteval.relimp(riBoot[[i]])@lmg.lower,
#     booteval.relimp(riBoot[[i]])@est[-1:-3]
#   )
#   colnames(relimpBootResult) <- gsub(".lmg", "", colnames(relimpBootResult))
#   rownames(relimpBootResult) <- c("Upper", "Lower", "M")
#   bp <- barplot(relimpBootResult['M',], xlab="", ylab="R² 貢獻量", beside=T, ylim=c(0,0.4), las=2)
#   segments(bp, relimpBootResult["Upper",], bp, relimpBootResult["Lower",], lwd=2)
#   # title(names(fit.o)[i])
# }
# # quartz.save( "slide/beta相對貢獻量boot.pdf", type="pdf")
# quartz.save( "slide/beta相對貢獻量.png", type="png", dpi=300)
# dev.off()



########### find best fitting BAF
dB0t.pred <- data.frame(
  PCA.1 = predict(fit.o$PCA.1, dG1t),
  PCA.F = predict(fit.o$PCA.F, dG1t),
  FA.1  = predict(fit.o$FA.1, dG1t),
  FA.F  = predict(fit.o$FA.F, dG1t),
  MDS.1 = predict(fit.o$MDS.1, dG1t),
  MDS.F = predict(fit.o$MDS.F, dG1t),
  RDA.1 = predict(fit.o$RDA.1, dG1t),
  RDA.F = predict(fit.o$RDA.F, dG1t)
)
quartz(width=9, height=9)
par(mfrow=c(1,4), cex=6/12)
pairs(
  cbind(dB0t.pred, dB0t),
  lower.panel=panel.smooth,
  upper.panel=panel.cor, pch=".")
quartz.save("slide/all.png", dpi=300, type="png")
# quartz.save("slide/all.pdf", type="pdf")
dev.off()



# 選擇表
rbind(
  `# IV` = sapply(fit.o, FUN = function(fit) {length(attr(fit$terms, "term.labels"))} ),
  NRMSE = sapply(fit.o, FUN = function(fit) {sqrt(sum(resid(fit)^2) / length(fit$residuals)) / sd(fit$model[[1]])} ),
  R2 = sapply(fit.o, FUN = function(fit) {summary(fit)$r.squared} ),
  `Adjusted R2` = sapply(fit.o, FUN = function(fit) {summary(fit)$adj.r.squared} ),
  `Kendall tau` = sapply(fit.o, FUN = function(fit) {
    a <- predict(fit, fit$model)
    b <- fit$model[[1]]
    return(cor(a,b,method="kendall"))
  }),
  `AIC/1000` = sapply(fit.o, AIC) / 1000,
  `BIC/1000` = sapply(fit.o, BIC) / 1000
)  %T>% print %>%
  xtable(., digits=3, align = c("l", rep("r", ncol(.))) ) %>%
  print(., NA.string="—", booktabs=F, size="small", math.style.negative=T,
    file="slide/選擇表.tex",
    sanitize.colnames.function = function(x){gsub("(PCA|FA|MDS|RDA).(AIC|BIC)", "\\1\\\\textsubscript{\\2}", x)  },
    sanitize.rownames.function = function(x){
        x[x=="(Intercept)"] <- "Int."
        x <- gsub("([A-Z])2$", "$ \\1^2 $", x)
        x <- sub("#", "\\\\#", x)
        return(x)
    },
    tabular.environment = "mytabular"
  )

quartz(width=7, height=2)
par(mar=c(2,4,0,0)+0.1, family = "Noto Sans CJK TC", mgp=c(3,1,0), cex=10/12)
boxplot(sapply(fit.o, resid) %>% scale(., center=F), xlab="模型", ylab="標準化殘差", las=1)
quartz.save("slide/resid最佳迴歸法.png", dpi=300, type="png")
# quartz.save("./slide/resid最佳迴歸法.pdf", type="pdf")
dev.off()
# quartz(width=7, height=2)
# tmp <- fit.o %>% sapply(., resid) %>% as.data.frame %>% gather(., Model, Residual, 1:8)
# ggplot(tmp, aes(Model, Residual)) + geom_jitter(colour=gray(0.5),size=1,height=10) + geom_violin(alpha=0.8) + theme_bw()




## beta plot before standardizing
# tmp <- data.frame(
#   aiccol = 1:4,
#   biccol = 5:8,
#   filename = paste0("beta-", sub(".AIC", "", colnames(coefM.o)[1:4]))
# )
# for(i in 1:nrow(tmp)) {
#   quartz(width=7, height=3.5)
#   par(mar=c(5,5,0,0)+0.1, family = "Noto Sans CJK TC", mgp=c(4,1,0), cex=8/12)
#   plot(NULL, xlim=c(1, nrow(coefM.o)), ylim=c(-15,15), xaxt="n", yaxt="n", xlab="類型", ylab="原始迴歸係數" )
#     axis(1, 1:nrow(coefM.o), rownames(coefM.o), las=2)
#     axis(2, seq(-16, 16, 2), seq(-16, 16, 2), las=2)
#     abline(h=0)
#     arrows(
#       1:nrow(coefM.o)-0.2,
#       rep(0, nrow(coefM.o)),
#       1:nrow(coefM.o)-0.2,
#       ifelse(
#         coefM.o[, tmp$aiccol[i]] >= 0 | is.na(coefM.o[, tmp$aiccol[i]]),
#         coefM.o[,tmp$aiccol[i]],
#         coefM.o[,tmp$aiccol[i]]
#       ),
#       code=2, length = 0.05, col=1 #, lwd=abs(coefM.o[,1])/max(abs(coefM.o[,1]))*5+1
#     )
#     arrows(
#       1:nrow(coefM.o)+0.2,
#       rep(0, nrow(coefM.o)),
#       1:nrow(coefM.o)+0.2,
#       ifelse(
#         coefM.o[,tmp$biccol[i]] >= 0 | is.na(coefM.o[,tmp$biccol[i]]),
#         coefM.o[,tmp$biccol[i]],
#         coefM.o[,tmp$biccol[i]]
#       ) ,
#       code=2, length = 0.05, col=2, lty=1#, lwd=abs(coefM.o[,2])/max(abs(coefM.o[,2]))*5+1
#     )
#     legend(1, 8, c("AIC", "BIC"), lty=c(1,1), col=c(1,2))
#   quartz.save( paste0(tmp$filename[i], ".png") , dpi=300, type="png")
#   quartz.save( paste0(tmp$filename[i], ".pdf"), type="pdf")
#   dev.off()
# }



















# 創造 10-fold 並回傳 trainning set
# 切出 TD & VD
set.seed(52004800)
# set.seed(1000)
# set.seed(124234545)
TD.row <- createFolds(1:nrow(dB0t), k = 10, list = TRUE, returnTrain = T)
# TD.row <- createMultiFolds(1:nrow(dB0t), k = 10, times=10)
VD.row <- lapply(TD.row, function(x){ which( !((1:nrow(dB0t)) %in% x) )})
dB0t.TD <- lapply(TD.row, function(x){dB0t[x,]})
dG1t.TD <- lapply(TD.row, function(x){dG1t[x,]})
dB0t.VD <- lapply(VD.row, function(x){dB0t[x,]})
dG1t.VD <- lapply(VD.row, function(x){dG1t[x,]})
dB0t.dG1t.TD <- apply(mapply(c, dB0t.TD, dG1t.TD), 2, function(this){as.data.frame(this)})
dB0t.dG1t.VD <- apply(mapply(c, dB0t.VD, dG1t.VD), 2, function(this){as.data.frame(this)})












#### 10-fold CV 找哪個模型 fit.o[[i]] 最好
#### 資料來自全資料，不再自己降維，再訓練出不同的迴歸式
#### 目標：求最佳迴歸方法（不是求最佳的已知迴歸式）
Bscore.TD <- list()
Bscore.VD <- list()
fit.TD <- list()
Bscore.TD.pred <- list()
for(i in 1:length(TD.row)){
  # 結合已降維的 dB0t.TD[[i]] 成 Bscore.TD[[i]]
  Bscore.TD[[i]] <- data.frame(NO = dt$NO[TD.row[[i]]])
  Bscore.TD[[i]]$PCA.1 <- Bscore.o[TD.row[[i]],]$PCA
  Bscore.TD[[i]]$PCA.F <- Bscore.o[TD.row[[i]],]$PCA
  Bscore.TD[[i]]$FA.1  <- Bscore.o[TD.row[[i]],]$FA
  Bscore.TD[[i]]$FA.F  <- Bscore.o[TD.row[[i]],]$FA
  Bscore.TD[[i]]$MDS.1 <- Bscore.o[TD.row[[i]],]$MDS
  Bscore.TD[[i]]$MDS.F <- Bscore.o[TD.row[[i]],]$MDS
  Bscore.TD[[i]]$RDA.1 <- Bscore.o[TD.row[[i]],]$RDA
  Bscore.TD[[i]]$RDA.F <- Bscore.o[TD.row[[i]],]$RDA
  rownames(Bscore.TD[[i]]) <- dt$NO[TD.row[[i]]]
  Bscore.TD[[i]]$NO <- NULL

  # 結合已降維的 dB0t.VD[[i]] 成 Bscore.VD[[i]]
  Bscore.VD[[i]] <- data.frame(NO = dt$NO[VD.row[[i]]])
  Bscore.VD[[i]]$PCA.1 <- Bscore.o[VD.row[[i]],]$PCA
  Bscore.VD[[i]]$PCA.F <- Bscore.o[VD.row[[i]],]$PCA
  Bscore.VD[[i]]$FA.1  <- Bscore.o[VD.row[[i]],]$FA
  Bscore.VD[[i]]$FA.F  <- Bscore.o[VD.row[[i]],]$FA
  Bscore.VD[[i]]$MDS.1 <- Bscore.o[VD.row[[i]],]$MDS
  Bscore.VD[[i]]$MDS.F <- Bscore.o[VD.row[[i]],]$MDS
  Bscore.VD[[i]]$RDA.1 <- Bscore.o[VD.row[[i]],]$RDA
  Bscore.VD[[i]]$RDA.F <- Bscore.o[VD.row[[i]],]$RDA

  rownames(Bscore.VD[[i]]) <- dt$NO[VD.row[[i]]]
  Bscore.VD[[i]]$NO <- NULL

  # multiple regression: Bscore.TD[[i]]$xxx ~ ., data=dG1t.TD
  fit.TD[[i]] <-  vector(mode = "list", 8)
  names(fit.TD[[i]]) <- paste( c("PCA", "PCA", "FA", "FA", "MDS", "MDS", "RDA", "RDA"), rep(c("1", "F"), 4), sep="." )
  for(j in 1:length(fit.TD[[i]])) {
    m1 <- lm(Bscore.TD[[i]] [[j]] ~ A+B+C+D+E+F+I+J , data=dG1t.TD[[i]])
    m0 <- lm(Bscore.TD[[i]] [[j]] ~ A+D+E+F+I , data=dG1t.TD[[i]])
    m2 <- lm(Bscore.TD[[i]] [[j]] ~ . , data=dG1t.TD[[i]])
    if(j %% 2 == 1) {
      fit.TD[[i]][[j]] <- stats::step(m1, scope = list(upper = formula(m2), lower = formula(m0)), k = this.k[j], direction = stepDirection, trace=1)
    } else {
      fit.TD[[i]][[j]] <- stats::step(m2, scope = list(upper = formula(m2), lower = formula(m0)), k = this.k[j], direction = stepDirection, trace=1)
    }
  }

  # 以訓練出的 fit.TD$xxx.yIC 模型，求 dG1t.VD 的預測值
  Bscore.TD.pred[[i]] <- data.frame(
    PCA.1 = predict(fit.TD[[i]]$PCA.1, dG1t.VD[[i]]),
    PCA.F = predict(fit.TD[[i]]$PCA.F, dG1t.VD[[i]]),
    FA.1  = predict(fit.TD[[i]]$FA.1, dG1t.VD[[i]]),
    FA.F  = predict(fit.TD[[i]]$FA.F, dG1t.VD[[i]]),
    MDS.1 = predict(fit.TD[[i]]$MDS.1, dG1t.VD[[i]]),
    MDS.F = predict(fit.TD[[i]]$MDS.F, dG1t.VD[[i]]),
    RDA.1 = predict(fit.TD[[i]]$RDA.1, dG1t.VD[[i]]),
    RDA.F = predict(fit.TD[[i]]$RDA.F, dG1t.VD[[i]])
  )
}
# 結合所有預測 Bscore.TD.pred.final 和所有驗證 Bscore.VD.final
Bscore.TD.pred.求最佳迴歸法  <- do.call("rbind", Bscore.TD.pred)
  Bscore.TD.pred.求最佳迴歸法 <- Bscore.TD.pred.求最佳迴歸法[order(as.numeric(rownames(Bscore.TD.pred.求最佳迴歸法))), ]
Bscore.VD.求最佳迴歸法 <- do.call("rbind", Bscore.VD)
  Bscore.VD.求最佳迴歸法 <- Bscore.VD.求最佳迴歸法[order(as.numeric(rownames(Bscore.VD.求最佳迴歸法))), ]
dim(Bscore.TD.pred.求最佳迴歸法) ; dim(Bscore.VD.求最佳迴歸法)
names(Bscore.TD.pred.求最佳迴歸法) ; names(Bscore.VD.求最佳迴歸法)
names(Bscore.TD.pred.求最佳迴歸法) <- sub(".p", "", names(Bscore.TD.pred.求最佳迴歸法))
# 看看誰最好
rbind(
  # NRMSE = sqrt(colSums((Bscore.TD.pred.求最佳迴歸法 - Bscore.VD.求最佳迴歸法)^2) / nrow(Bscore.TD.pred.求最佳迴歸法)) / apply(Bscore.VD.求最佳迴歸法, 2, sd)
  NRMSE = rmse(Bscore.TD.pred.求最佳迴歸法, Bscore.VD.求最佳迴歸法) / apply(Bscore.VD.求最佳迴歸法, 2, sd),
  `Pearson r` = cor(Bscore.TD.pred.求最佳迴歸法, Bscore.VD.求最佳迴歸法) %>% diag,
  `Kendall tau` = cor(Bscore.TD.pred.求最佳迴歸法, Bscore.VD.求最佳迴歸法, method="kendall") %>% diag
  )  %T>% print %>%
    xtable(., digits=3, align = c("l", rep("r", ncol(.))) ) %>%
    print(., NA.string="—", booktabs=F, size="small", math.style.negative=T,
      file="slide/CV求最佳迴歸法.tex",
      sanitize.colnames.function = function(x){gsub("(PCA|FA|MDS|RDA).(1|F)", "\\1\\\\textsubscript{\\2}", x)  },
      sanitize.rownames.function = function(x){
          x[x=="(Intercept)"] <- "Int."
          x <- gsub("([A-Z])2$", "$ \\1^2 $", x)
          x <- sub("#", "\\\\#", x)
          return(x)
      },
      tabular.environment = "mytabular"
    )
quartz(width=7, height=2)
par(mar=c(2,4,0,0)+0.1, family = "Noto Sans CJK TC", mgp=c(3,1,0), cex=10/12)
boxplot(  Bscore.TD.pred.求最佳迴歸法 - Bscore.VD.求最佳迴歸法 , xlab="模型", ylab="訓練與驗證差值", las=1)
quartz.save("slide/CV求最佳迴歸法.png", dpi=300, type="png")
# quartz.save("slide/CV求最佳迴歸法.pdf", type="pdf")
dev.off()




















#### 10-fold CV 找8種迴歸因子 哪組變數組合最好
#### 不是挑迴歸係數，而是挑因子
formu.o <- list(
  PCA.AIC = as.character((as.list(fit.o[[1]]$call)$formula)[[3]])[2],
  FA.AIC  = as.character((as.list(fit.o[[2]]$call)$formula)[[3]])[2],
  MDS.AIC = as.character((as.list(fit.o[[3]]$call)$formula)[[3]])[2],
  RDA.AIC = as.character((as.list(fit.o[[4]]$call)$formula)[[3]])[2],
  PCA.BIC = as.character((as.list(fit.o[[5]]$call)$formula)[[3]])[2],
  FA.BIC  = as.character((as.list(fit.o[[6]]$call)$formula)[[3]])[2],
  MDS.BIC = as.character((as.list(fit.o[[7]]$call)$formula)[[3]])[2],
  RDA.BIC = as.character((as.list(fit.o[[8]]$call)$formula)[[3]])[2]
)
Bscore.TD <- list()
Bscore.VD <- list()
fit.TD <- list()
Bscore.TD.pred <- list()
for(i in 1:length(TD.row)){
  # 生態訓練資料 dB0t.TD[[i]] 重新降維成 Bscore.TD[[i]] (x)
  # 已降維生態訓練資料 dB0t.TD[[i]] 命為 Bscore.TD[[i]] (o)
  Bscore.TD[[i]] <- data.frame(NO = dt$NO[TD.row[[i]]])
  Bscore.TD[[i]]$PCA.1 <- Bscore.o[TD.row[[i]],]$PCA
  Bscore.TD[[i]]$PCA.F <- Bscore.o[TD.row[[i]],]$PCA
  Bscore.TD[[i]]$FA.1  <- Bscore.o[TD.row[[i]],]$FA
  Bscore.TD[[i]]$FA.F  <- Bscore.o[TD.row[[i]],]$FA
  Bscore.TD[[i]]$MDS.1 <- Bscore.o[TD.row[[i]],]$MDS
  Bscore.TD[[i]]$MDS.F <- Bscore.o[TD.row[[i]],]$MDS
  Bscore.TD[[i]]$RDA.1 <- Bscore.o[TD.row[[i]],]$RDA
  Bscore.TD[[i]]$RDA.F <- Bscore.o[TD.row[[i]],]$RDA
  rownames(Bscore.TD[[i]]) <- dt$NO[TD.row[[i]]]
  Bscore.TD[[i]]$NO <- NULL

  # 生態驗證資料 dB0t.VD[[i]] = Bscore.VD[[i]]          (x)
  # 已降維生態驗證資料 dB0t.VD[[i]] 命為 Bscore.VD[[i]] (o)
  Bscore.VD[[i]] <- data.frame(NO = dt$NO[VD.row[[i]]])
  Bscore.VD[[i]]$PCA.1 <- Bscore.o[VD.row[[i]],]$PCA
  Bscore.VD[[i]]$PCA.F <- Bscore.o[VD.row[[i]],]$PCA
  Bscore.VD[[i]]$FA.1  <- Bscore.o[VD.row[[i]],]$FA
  Bscore.VD[[i]]$FA.F  <- Bscore.o[VD.row[[i]],]$FA
  Bscore.VD[[i]]$MDS.1 <- Bscore.o[VD.row[[i]],]$MDS
  Bscore.VD[[i]]$MDS.F <- Bscore.o[VD.row[[i]],]$MDS
  Bscore.VD[[i]]$RDA.1 <- Bscore.o[VD.row[[i]],]$RDA
  Bscore.VD[[i]]$RDA.F <- Bscore.o[VD.row[[i]],]$RDA
  rownames(Bscore.VD[[i]]) <- dt$NO[VD.row[[i]]]
  Bscore.VD[[i]]$NO <- NULL

  # 拿 Bscore.TD[[i]]$xxA.XIC 於 formu.o$xxA.xIC 做複迴歸
  # 重點：迴歸式的IV組合已預先確定
  fit.TD[[i]] <-  vector(mode = "list", 8)
  names(fit.TD[[i]]) <- paste( c("PCA", "PCA", "FA", "FA", "MDS", "MDS", "RDA", "RDA"), rep(c("F", "1"), 4), sep="." )
  for(j in 1:length(fit.TD[[i]])) {
    fit.TD[[i]][[j]] <- lm(
      as.formula(paste(c( "Bscore.TD[[i]][[j]] ~ ", formu.o[[j]] ), collapse="")),
      data=dG1t.TD[[i]]
    )
  }

  # 以 fit.TD[[i]]$xxx.yIC 模型，求 dG1t.TD 的預測值 Bscore.TD.pred
  Bscore.TD.pred[[i]] <- data.frame(
    PCA.1 = predict(fit.TD[[i]]$PCA.1, dG1t.VD[[i]]),
    PCA.F = predict(fit.TD[[i]]$PCA.F, dG1t.VD[[i]]),
    FA.1  = predict(fit.TD[[i]]$FA.1, dG1t.VD[[i]]),
    FA.F  = predict(fit.TD[[i]]$FA.F, dG1t.VD[[i]]),
    MDS.1 = predict(fit.TD[[i]]$MDS.1, dG1t.VD[[i]]),
    MDS.F = predict(fit.TD[[i]]$MDS.F, dG1t.VD[[i]]),
    RDA.1 = predict(fit.TD[[i]]$RDA.1, dG1t.VD[[i]]),
    RDA.F = predict(fit.TD[[i]]$RDA.F, dG1t.VD[[i]])
  )
}
# 結合所有預測 Bscore.TD.pred.final 和所有驗證 Bscore.VD.final
Bscore.TD.pred.最佳解釋變數組 <- do.call("rbind", Bscore.TD.pred)
  Bscore.TD.pred.最佳解釋變數組 <- Bscore.TD.pred.最佳解釋變數組[order(as.numeric(rownames(Bscore.TD.pred.最佳解釋變數組))), ]
Bscore.VD.最佳解釋變數組 <- do.call("rbind", Bscore.VD)
  Bscore.VD.最佳解釋變數組 <- Bscore.VD.最佳解釋變數組[order(as.numeric(rownames(Bscore.VD.最佳解釋變數組))), ]
dim(Bscore.TD.pred.最佳解釋變數組) ; dim(Bscore.VD.最佳解釋變數組)
names(Bscore.TD.pred.最佳解釋變數組) ; names(Bscore.VD.最佳解釋變數組)
names(Bscore.TD.pred.最佳解釋變數組) <- sub(".p", "", names(Bscore.TD.pred.最佳解釋變數組))
# 看看誰最好
# rbind(NRMSE = sqrt(colSums((Bscore.TD.pred.最佳解釋變數組 - Bscore.VD.最佳解釋變數組)^2) / nrow(Bscore.TD.pred.最佳解釋變數組)) / apply(Bscore.VD.最佳解釋變數組, 2, sd))  %>% print %>% xtable(., digits=4) %>% print.xtable(., booktabs=T)
rbind(
  NRMSE = rmse(Bscore.TD.pred.最佳解釋變數組, Bscore.VD.最佳解釋變數組) / apply(Bscore.VD.最佳解釋變數組, 2, sd),
  `Pearson r` = cor(Bscore.TD.pred.最佳解釋變數組, Bscore.VD.最佳解釋變數組) %>% diag,
  `Kendall tau` = cor(Bscore.TD.pred.最佳解釋變數組, Bscore.VD.最佳解釋變數組, method="kendall") %>% diag
  )  %T>% print %>%
    xtable(., digits=3, align = c("l", rep("r", ncol(.))) ) %>%
    print(., NA.string="—", booktabs=F, size="small", math.style.negative=T,
      file="slide/CV-最佳解釋變數組.tex",
      sanitize.colnames.function = function(x){gsub("(PCA|FA|MDS|RDA).(1|F)", "\\1\\\\textsubscript{\\2}", x)  },
      sanitize.rownames.function = function(x){
          x[x=="(Intercept)"] <- "Int."
          x <- gsub("([A-Z])2$", "$ \\1^2 $", x)
          x <- sub("#", "\\\\#", x)
          return(x)
      },
      tabular.environment = "mytabular"
    )
quartz(width=7, height=2)
par(mar=c(2,4,0,0)+0.1, family = "Noto Sans CJK TC", mgp=c(3,1,0), cex=10/12)
boxplot(  Bscore.TD.pred.最佳解釋變數組 - Bscore.VD.最佳解釋變數組 , xlab="模型", ylab="訓練與驗證差值", lax=1)
quartz.save("slide/CV-最佳解釋變數組.png", dpi=300, type="png")
# quartz.save("slide/CV-最佳解釋變數組.pdf", type="pdf")
dev.off()











#### 10-fold CV 找8種求解法何者較佳
#### 目標：求最佳整套求解方式（包括降維與迴歸方式）
mds.try.n <- 20
Bscore.TD <- list()
Bscore.VD <- list()
fit.TD <- list()
Bscore.TD.pred <- list()
for(i in 1:length(TD.row)){
  # 生態訓練資料 dB0t.TD[[i]] 重新降維成 Bscore.TD[[i]]
  Bscore.TD[[i]] <- data.frame(NO = dt$NO[TD.row[[i]]])
  Bscore.TD[[i]]$PCA.1 <- fitPCA(dB0t.TD[[i]])
  Bscore.TD[[i]]$PCA.F <- Bscore.TD[[i]]$PCA.1
  Bscore.TD[[i]]$FA.1  <- fitFA(dB0t.TD[[i]])
  Bscore.TD[[i]]$FA.F  <- Bscore.TD[[i]]$FA.1
  Bscore.TD[[i]]$MDS.1 <- fitMDS(dB0t.TD[[i]])
  Bscore.TD[[i]]$MDS.F <- Bscore.TD[[i]]$MDS.1
  Bscore.TD[[i]]$RDA.1 <- fitRDA(dB0t.TD[[i]], dG1t.TD[[i]])
  Bscore.TD[[i]]$RDA.F <- Bscore.TD[[i]]$RDA.1
  rownames(Bscore.TD[[i]]) <- dt$NO[TD.row[[i]]]
  Bscore.TD[[i]]$NO <- NULL

  # 生態驗證資料 dB0t.VD[[i]] = Bscore.VD[[i]]
  Bscore.VD[[i]] <- data.frame(NO = dt$NO[VD.row[[i]]])
  Bscore.VD[[i]]$PCA.1 <- fitPCA(dB0t.VD[[i]])
  Bscore.VD[[i]]$PCA.F <- Bscore.VD[[i]]$PCA.1
  Bscore.VD[[i]]$FA.1  <- fitFA(dB0t.VD[[i]])
  Bscore.VD[[i]]$FA.F  <- Bscore.VD[[i]]$FA.1
  Bscore.VD[[i]]$MDS.1 <- fitMDS(dB0t.VD[[i]])
  Bscore.VD[[i]]$MDS.F <- Bscore.VD[[i]]$MDS.1
  Bscore.VD[[i]]$RDA.1 <- fitRDA(dB0t.VD[[i]], dG1t.VD[[i]])
  Bscore.VD[[i]]$RDA.F <- Bscore.VD[[i]]$RDA.1
  rownames(Bscore.VD[[i]]) <- dt$NO[VD.row[[i]]]
  Bscore.VD[[i]]$NO <- NULL

  # multiple regression: Bscore.TD[[i]]$xxx ~ ., data=dG1t.TD
  fit.TD[[i]] <-  vector(mode = "list", 8)
  names(fit.TD[[i]]) <- paste( c("PCA", "PCA", "FA", "FA", "MDS", "MDS", "RDA", "RDA"), rep(c("1", "F"), 4), sep="." )
  this.k <- rep(c(2, 2), 4)
  for(j in 1:length(fit.TD[[i]])) {
    m1 <- lm(Bscore.TD[[i]] [[j]] ~ A+B+C+D+E+F+I+J , data=dG1t.TD[[i]])
    m0 <- lm(Bscore.TD[[i]] [[j]] ~ 1 , data=dG1t.TD[[i]])
    m2 <- lm(Bscore.TD[[i]] [[j]] ~ . , data=dG1t.TD[[i]])
    if(j %% 2 == 1) {
      # 從 reduced model m1 開始, AIC
      fit.TD[[i]][[j]] <- stats::step(m1, scope = list(upper = formula(m2), lower = formula(m0)), k = this.k[j], direction = stepDirection, trace=1)
    } else {
      # 從 full model m2 開始, BIC
      fit.TD[[i]][[j]] <- stats::step(m2, scope = list(upper = formula(m2), lower = formula(m0)), k = this.k[j], direction = stepDirection, trace=1)
    }
  }

  # 以 fit.TD[[i]]$xxx.yIC 模型，求 dG1t.TD 的預測值 Bscore.TD.pred
  Bscore.TD.pred[[i]] <- data.frame(
    PCA.1 = predict(fit.TD[[i]]$PCA.1, dG1t.VD[[i]]),
    PCA.F = predict(fit.TD[[i]]$PCA.F, dG1t.VD[[i]]),
    FA.1  = predict(fit.TD[[i]]$FA.1, dG1t.VD[[i]]),
    FA.F  = predict(fit.TD[[i]]$FA.F, dG1t.VD[[i]]),
    MDS.1 = predict(fit.TD[[i]]$MDS.1, dG1t.VD[[i]]),
    MDS.F = predict(fit.TD[[i]]$MDS.F, dG1t.VD[[i]]),
    RDA.1 = predict(fit.TD[[i]]$RDA.1, dG1t.VD[[i]]),
    RDA.F = predict(fit.TD[[i]]$RDA.F, dG1t.VD[[i]])
  )
}
# 結合所有預測 Bscore.TD.pred.final 和所有驗證 Bscore.VD.final
# NO==640 八卦山7很怪
Bscore.TD.pred.final <- do.call("rbind", Bscore.TD.pred)
  Bscore.TD.pred.final <- Bscore.TD.pred.final[order(as.numeric(rownames(Bscore.TD.pred.final))), ]
Bscore.VD.final <- do.call("rbind", Bscore.VD)
  Bscore.VD.final <- Bscore.VD.final[order(as.numeric(rownames(Bscore.VD.final))), ]
dim(Bscore.TD.pred.final) ; dim(Bscore.VD.final)
names(Bscore.TD.pred.final) ; names(Bscore.VD.final)
names(Bscore.TD.pred.final) <- sub(".p", "", names(Bscore.TD.pred.final))
# 看看誰最好
rbind(
  NRMSE = rmse(Bscore.TD.pred.final, Bscore.VD.final) / apply(Bscore.VD.final, 2, sd),
  `Pearson r` = cor(Bscore.TD.pred.final, Bscore.VD.final) %>% diag,
  `Kendall tau` = cor(Bscore.TD.pred.final, Bscore.VD.final, method="kendall") %>% diag
  )  %T>% print %>%
    xtable(., digits=3, align = c("l", rep("r", ncol(.))) ) %>%
    print(., NA.string="—", booktabs=F, size="small", math.style.negative=T,
      file="slide/CV-final.tex",
      sanitize.colnames.function = function(x){gsub("(PCA|FA|MDS|RDA).(1|F)", "\\1\\\\textsubscript{\\2}", x)  },
      sanitize.rownames.function = function(x){
          x[x=="(Intercept)"] <- "Int."
          x <- gsub("([A-Z])2$", "$ \\1^2 $", x)
          x <- sub("#", "\\\\#", x)
          return(x)
      },
      tabular.environment = "mytabular"
    )
quartz(width=7, height=2)
par(mar=c(2,4,0,0)+0.1, family = "Noto Sans CJK TC", mgp=c(3,1,0), cex=10/12)
boxplot(  Bscore.TD.pred.final - Bscore.VD.final , xlab="模型", ylab="訓練與驗證差值", lax=1)
quartz.save("slide/CV-final.png", dpi=300, type="png")
# quartz.save("slide/CV-final.pdf", type="pdf")
dev.off()
# max(Bscore.TD.pred.final[,1])
# dt [max(Bscore.TD.pred.final[,1]) == Bscore.TD.pred.final[,1] , ]
# kruskal.test(value ~ Var2,  data=melt(tmp))
# posthoc.kruskal.nemenyi.test(melt(tmp)$value, melt(tmp)$Var2, method="Tukey")




########## 印出迴歸式
coefM.o[c("FA.1", "RDA.1")] %>% na.omit %>%
  xtable(., digits=5, align = c("l", rep("r", ncol(.))) ) %>%
  print(., NA.string="", booktabs=F, size="normalsize", math.style.negative=T,
    file="slide/印出迴歸式.tex",
    sanitize.colnames.function = function(x){gsub("(PCA|FA|MDS|RDA).(AIC|BIC)", "\\1\\\\textsubscript{\\2}", x)  },
    sanitize.rownames.function = function(x){
        x[x=="(Intercept)"] <- "Int."
        x <- gsub("^([A-Z])$", "$ \\1 $", x)
        x <- gsub("^([A-Z])2$", "$ \\1^2 $", x)
        x <- gsub("^([A-Z])([A-Z])$", " $\\1 \\\\times \\2 $", x)
        return(x)
    },
    tabular.environment = "mytabular"
  )





############ 檢視預測
預測排名.FA.1 <- data.frame(
  dt[, colnames(dt) %in% c("Location", "Spot", LETTERS[1:26])],
  Rank木 = rank(-dt$HsWood),
  Rank草 = rank(-dt$HsGrass),
  Rank蛛 = rank(-dt$HsSpider),
  Rank蟲 = rank(-dt$HsInsect),
  預測名次 = rank(-dB0t.pred$FA.1)
)  %>%
  .[, !colnames(.) %in% c("G", "H", "K", "L", "NO", "DgsWood", "DgsGrass", "DgsSpider", "DgsInsect")]
預測排名.FA.1 %>%
  .[order(.$預測名次) , ] %>%
  .[1:30 , ] %>%
  xtable(., digits = c(rep(0,13), rep(0,4), 0) ) %>%
  print(., tabular.environment = "mytabular", size="tiny", include.rownames=F, file="slide/預測1.tex")
預測排名.FA.1 %>%
  .[order(-.$預測名次) , ] %>%
  .[1:30 , ] %>%
  xtable(., digits = c(rep(0,13), rep(0,4), 0) ) %>%
  print(., tabular.environment = "mytabular", size="tiny", include.rownames=F, file="slide/預測2.tex")

預測排名.RDA.1 <- data.frame(
  dt[, colnames(dt) %in% c("Location", "Spot", LETTERS[1:26])],
  Rank木 = rank(-dt$HsWood),
  Rank草 = rank(-dt$HsGrass),
  Rank蛛 = rank(-dt$HsSpider),
  Rank蟲 = rank(-dt$HsInsect),
  預測名次 = rank(-dB0t.pred$RDA.1)
)  %>%
  .[, !colnames(.) %in% c("G", "H", "K", "L", "NO", "DgsWood", "DgsGrass", "DgsSpider", "DgsInsect")]
預測排名.RDA.1 %>%
  .[order(.$預測名次) , ] %>%
  .[1:30 , ] %>%
  xtable(., digits = c(rep(0,13), rep(0,4), 0) ) %>%
  print(., tabular.environment = "mytabular", size="tiny", include.rownames=F, file="slide/預測3.tex")
預測排名.RDA.1 %>%
  .[order(-.$預測名次) , ] %>%
  .[1:30 , ] %>%
  xtable(., digits = c(rep(0,13), rep(0,4), 0) ) %>%
  print(., tabular.environment = "mytabular", size="tiny", include.rownames=F, file="slide/預測4.tex")




#### 10-fold CV 找哪種降維方法最好 ----
#### 不可能，因為無法訓練出 FA MDS 的模型來預測分數
#### 改做 bootstrap kendall correlation 找哪種降維方法最好
# k <- 10000
# ind <- vector("list", k)
# set.seed(52004800)#; set.seed(19821006); set.seed(124234545)
# ind <- lapply(ind, function(x){return(sample(1:nrow(dB0t), replace=T))})
# bootR <- vector("list", k)
# pb <- txtProgressBar(min = 0, max = k, style = 3)
# for(i in 1:k){
#   bootR[[i]] <- data.frame(
#     cor(Bscore.o[ind[[i]],-1]$PCA, dB0t[ind[[i]],], method="kendall" ),
#     cor(Bscore.o[ind[[i]],-1]$FA,  dB0t[ind[[i]],], method="kendall" ),
#     cor(Bscore.o[ind[[i]],-1]$MDS, dB0t[ind[[i]],], method="kendall" ),
#     cor(Bscore.o[ind[[i]],-1]$RDA, dB0t[ind[[i]],], method="kendall" )
#   )
#   setTxtProgressBar(pb, i)
# }
# close(pb)
# bootR.final <- do.call("rbind", bootR)
# colnames(bootR.final) <- paste0(
#   rep(c("HsWood", "HsGrass", "HsSpider", "HsInsect"), 4),
#   c(rep(".PCA", 4), rep(".FA", 4), rep(".MDS", 4), rep(".RDA", 4))
# )
# apply(bootR.final, 2, function(x){
#     quantile(x, c(0.0005,0.005,0.025,0.975,0.995,0.9995))
# })
# save(bootR.final, file="boorT.final.Rdata")










######## 回推未中心化的迴歸式
round(coefM.o[c("PCA.AIC","RDA.AIC")], digits=10)
finalFormula <- vector("list", 2)
names(finalFormula) <- c("PCA.AIC","RDA.AIC")

# dG1t1Mean
i = 1
finalFormula[[i]] <- list()
finalFormula[[i]]$intercept <-
  coefM.o["(Intercept)", i] +
  -coefM.o["A", i] * dG1t1Mean["A"] +
  -coefM.o["B", i] * dG1t1Mean["B"] +
  -coefM.o["E", i] * dG1t1Mean["E"] +
  -coefM.o["J", i] * dG1t1Mean["J"] +
  coefM.o["A2", i] * dG1t1Mean["A"]^2 +
  coefM.o["AD", i] * dG1t1Mean["A"] * dG1t1Mean["D"] +
  coefM.o["AI", i] * dG1t1Mean["A"] * dG1t1Mean["I"] +
  coefM.o["BE", i] * dG1t1Mean["B"] * dG1t1Mean["E"] +
  coefM.o["BI", i] * dG1t1Mean["B"] * dG1t1Mean["I"] +
  coefM.o["DI", i] * dG1t1Mean["D"] * dG1t1Mean["I"]
finalFormula[[i]]$A <-
  coefM.o["A", i] +
  -2 * coefM.o["A2", i] * dG1t1Mean["A"] +
  -coefM.o["AD", i] * dG1t1Mean["D"] +
  -coefM.o["AI", i] * dG1t1Mean["I"]
finalFormula[[i]]$A2 <-
  coefM.o["A2", i]
finalFormula[[i]]$B <-
  coefM.o["B", i] +
  -coefM.o["BE", i] * dG1t1Mean["E"] +
  -coefM.o["BI", i] * dG1t1Mean["I"]
finalFormula[[i]]$D <-
  -coefM.o["AD", i] * dG1t1Mean["A"] +
  -coefM.o["DI", i] * dG1t1Mean["I"]
finalFormula[[i]]$E <-
  coefM.o["E", i] +
  -coefM.o["BE", i] * dG1t1Mean["B"]
finalFormula[[i]]$I <-
  -coefM.o["AI", i] * dG1t1Mean["A"] +
  -coefM.o["BI", i] * dG1t1Mean["B"] +
  -coefM.o["DI", i] * dG1t1Mean["D"]
finalFormula[[i]]$J <-
  coefM.o["J", i] +
  -coefM.o["AI", i] * dG1t1Mean["A"]
finalFormula.PCA.AIC <- do.call("c", finalFormula[[1]])
t(finalFormula.PCA.AIC) %>% xtable

  coefM.o["A2", 1] * dG1t1Mean["A"] +



#
#
# ############# 做出所有可能假資料
# areaL <- 4 # 分 areaL 級面積等級
# typeN <- 7  # 有 typeN 種面積
# sizeM <- areaL ^ (typeN)
# tmp0 <- vector("list", typeN)
# names(tmp0) <- c("A","B","C","D","E","I","J")
# for(i in 1:typeN){
#   tmp0[[i]] <- as.numeric(gl(areaL, areaL^(i-1), sizeM)) - 1
# }
# dGsimAll1 <- do.call("cbind", tmp0) %>%
#   data.frame(., check.names = F, stringsAsFactors = F) %>%
#   decostand(., "total", 1) %>%
#   unique %>%
#   apply(., 2, function(x){x-mean(x)}) %>%
#   as.data.frame
# # dGsimAll1 完成：已建好一次式
# dGsimAll2 <- as.data.frame(dGsimAll1^2); names(dGsimAll2) <- paste0(names(dGsimAll2), "2")
# # dGsimAll2 完成：已建好二次式和交互作用
# dGsimAll <- cbind(as.data.frame(model.matrix(~(.)^2, dGsimAll1)), dGsimAll2)  [,-1]
# names(dGsimAll) <- gsub(":", "", names(dGsimAll))
# dGsimAll <- dGsimAll[
#   c(names(fit.o$PCA.AIC$model)[-1], names(fit.o$RDA.AIC$model)[-1]) %>% unique
# ] # 只挑出 fit.o$RDA.AIC 用得到的IV變數
# dGsimAll <- dGsimAll[, order(names(dGsimAll))]
# head(dGsimAll);dim(dGsimAll)
# # dGsimAll 完成：結合一次式、二次式和交互作用
# #### 代入迴歸式
# pred.dGsimAll.pca.aic <- predict(fit.o$PCA.AIC, dGsimAll)
# pred.dGsimAll.rda.aic <- predict(fit.o$RDA.AIC, dGsimAll)
#
# ## heatmap PCA.AIC
# x <- data.frame(
#   PCA.AIC = pred.dGsimAll.pca.aic,
#   RDA.AIC = pred.dGsimAll.rda.aic,
#   dGsimAll) %>%
#   as.matrix %>%
#   .[order(pred.dGsimAll.pca.aic),]
# quartz(width=7, height=7)
# par(family = "Noto Sans CJK TC", cex=10/12)
# hv <- heatmap(
#   #x,  # 很難用
#   x %>% apply(., 2, function(x){rank(x)}), #把 x rank 化比較好畫
#   Rowv = NA, Colv = NA,
#   # col = cm.colors(50),
#   col = rainbow(200, start=1/256, end=200/256),
#   # col = gray(seq(1, 0, length=10)),
#   # scale = "column",
#   scale = "none",
#   na.rm = F,
#   margins = c(5,22), cexRow=1, cexCol=1,
#   labRow = "",
#   # labCol = c("未定義", LETTERS[1:13]),
#   xlab = "", ylab =  "")
# # quartz.save("./slide/所有情況模擬PcaAic.pdf", type="pdf")
# quartz.save("slide/所有情況模擬PcaAic.png", dpi=300, type="png")
# dev.off()
#
#
#
#














# ##################### 亂數模型
# #### 創造假資料
# rown <- 500000    # 5000列假點
# coln <- 13      # 13類假面積
# set.seed(52004800)
# # 一次項
# dG1t1.Fake <- runif(rown * coln) %>%
#   matrix(., ncol = coln, nrow = rown) %>%
#   decostand(., "total", 1) %>%
#   as.data.frame
#   names(dG1t1.Fake) <- LETTERS[1:13]
# dG1t1.Fake <- dG1t1.Fake - matrix(rep(dG1t1Mean, rown), ncol=coln, byrow=T)
# # 二次項
# dG1t2.Fake <- dG1t1.Fake^2 %>% as.data.frame
#   names(dG1t2.Fake) %<>% paste0(., "2")
# # 交互作用項
# dG1tI.Fake <- model.matrix(~(.)^2, dG1t1.Fake) %>%
#   as.data.frame %>%
#   .[, grep(":", colnames(.))] %>%
#   as.data.frame
#   names(dG1tI.Fake) %<>% gsub(":", "", .)
# # 結合
# dG1t.Fake <- cbind(dG1t1.Fake, dG1t2.Fake, dG1tI.Fake) %>%
#   as.data.frame %>%
#   .[, colSums(.) != 0] #%>%   .[!colnames(.) %in% c("C2","F2","AC","AF","BC","BF","CD","CE","CF","CI","CJ","DF","EF","EI","EJ","FI","FJ")]
# str(dG1t.Fake)
# #### 代入迴歸式
# pred.pca.aic <- predict(fit.o$PCA.AIC, dG1t.Fake)
# pred.rda.aic <- predict(fit.o$RDA.AIC, dG1t.Fake)
# #### Top 20 sites pca.aic
# t.first <- cbind(
#   Bio = pred.pca.aic[order(-pred.pca.aic)[1:20]],
#   dG1t.Fake[order(-pred.pca.aic)[1:20], ] %>%
#     .[colnames(.) %in% LETTERS[1:13]] + matrix(rep(dG1t1Mean, rown), ncol=coln, byrow=T)
# )
# # #### Top -20 sites
# t.last <- cbind(
#   Bio = pred.pca.aic[order(pred.pca.aic)[1:20]],
#   dG1t.Fake[order(pred.pca.aic)[1:20], ] %>%
#     .[colnames(.) %in% LETTERS[1:13]] + matrix(rep(dG1t1Mean, rown), ncol=coln, byrow=T)
# )
# xtable(round(t.first, 2))
# xtable(round(t.last,  2))
