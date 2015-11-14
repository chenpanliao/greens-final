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
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
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

library(ggplot2)
quartzFonts(sans = quartzFont(rep("Noto Sans CJK TC Regular", 4)),
            serif = quartzFont(rep("Noto Sans CJK TC Regular", 4)))
BAF.mat <- matrix(
  c(1.4, 4.4, 3.9, 4.9, 4.3, 4.5, 3.7, 6.1, 7, 6.4, 4.3, 3.4, 5.4), 1,
  dimnames = list(c("BAF"), letters[1:13]) )
d <- read.csv("All_Data_20151020.csv")
names(d)


# 所有樣本
# 有 NA
dG0 <- d[, c(4:17)]
dG1 <- d[, c(5:17)]
dB0 <- d[, 18:25]


## 找出總面積不正常的點
quartz(width=7, height=3)
par(cex=10/12, mar=c(4,4,0,0)+0.1, family = "Noto Sans CJK TC")
hist(rowSums(dG0), nclass=40, main="", xlim=c(0,60000), ylim=c(0,500), col=3, xlab="總面積（m²）", ylab="次數")
  abline(v=c(8000,12000), lty=2, col=2)
  text(4000, 260, "< 8000\n(219/1169)")
  text(16000, 260, "> 12000\n(26/1169)")
quartz.save("invalid-area.pdf", type="pdf", family = "Noto Sans CJK TC")
quartz.save("invalid-area.png", dpi=300, type="png", family = "Noto Sans CJK TC")
# stem(rowSums(dG0), scale=2, width=100)
# table(rowSums(dG0) < 8000 | rowSums(dG0) > 12000)
dev.off()


## 找出未定義面積不正常的點
quartz(width=7, height=3)
par(cex=10/12, mar=c(4,4,0,0)+0.1, family = "Noto Sans CJK TC")
hist(dG0$Gx, nclass=40, main="", xlim=c(0,20000), ylim=c(0,600), col=3, xlab="未定義面積（m²）", ylab="次數")
  abline(v=c(8000), lty=2, col=2)
  # text(4000, 260, "< 8000\n(18.7%)")
  text(12000, 550, "> 8000\n(125/1169)")
quartz.save("invalid-gx-area.pdf", type="pdf", family = "Noto Sans CJK TC")
quartz.save("invalid-gx-area.png", dpi=300, type="png", family = "Noto Sans CJK TC")
# stem(rowSums(dG0), scale=2, width=100)
# table(rowSums(dG0) < 8000 | rowSums(dG0) > 12000)
dev.off()


## 有問題的景觀資料
x  <- as.matrix(d[, 4:17])
Location <- as.character(d$Location)
Location[duplicated(Location)] <- ""
quartz(width=7, height=7)
par(family = "Noto Sans CJK TC", cex=10/12)
hv <- heatmap(x,
  Rowv = NA, Colv = NA,
  # col = terrain.colors(50),
  col = gray(seq(1, 0, -0.01)),
  scale = "row", na.rm = F,
  margins = c(6,6),cexRow=0.3, cexCol=1,
  labRow = Location,
  labCol = c("未定義", LETTERS[1:13]),
  xlab = "類型", ylab =  "樣點")
quartz.save("invalid-type.pdf", type="pdf", family = "Noto Sans CJK TC")
quartz.save("invalid-type.png", dpi=300, type="png", family = "Noto Sans CJK TC")
dev.off()


## 有問題的多樣性資料
x  <- as.matrix(dB0)
Location <- as.character(d$Location)
Location[duplicated(Location)] <- ""
quartz(width=7, height=7)
par(family = "Noto Sans CJK TC", cex=10/12)
hv <- heatmap(x,
  Rowv = NA, Colv = NA,
  col = gray(seq(1, 0, -0.01)),
  scale = "column", na.rm = F,
  margins = c(6,6),cexRow=0.3,cexCol=1,
  labRow = Location,
  labCol = paste(c("木本","草本","蜘蛛","昆蟲"), c(rep("Simpson",4), rep("Shannon",4)), sep="\n"),
  xlab = "生物多樣性指數", ylab =  "樣點")
quartz.save("invalid-bio.pdf", type="pdf", family = "Noto Sans CJK TC")
quartz.save("invalid-bio.png", dpi=300, type="png", family = "Noto Sans CJK TC")
dev.off()















# 求解用的 dataset
# 無 NA
# 多樣性資料只看 shannon
# 去除 Gl
# 去除 < 7000 | > 13000 area
# 去除 Gx > 7000
# 面積以列標準化為 sum = 1
# NO=="640" 很怪
library(vegan)
dt <- d
dt <- subset(d, NO!="640")
dt <- na.omit(dt); dim(dt)
dt <- subset(dt, Gx+Ga+Gb+Gc+Gd+Ge+Gf+Gg+Gh+Gi+Gj+Gk+Gl+Gm >= 7000) ; dim(dt)
dt <- subset(dt, Gx+Ga+Gb+Gc+Gd+Ge+Gf+Gg+Gh+Gi+Gj+Gk+Gl+Gm <= 13000) ; dim(dt)
dt <- subset(dt, Gx <= 7000) ; dim(dt)
# 環境資料
dG0t1 <- dt[, c(4:17)]
dG0t2 <- dG0t1[, c(-8, -12,-13,-14)]
names(dG0t2)
dG0t3 <- decostand(dG0t2, "total", 1)
names(dG0t3);rowSums(dG0t3)
dG1t4 <- dG0t3[, c(-1, -8, -12,-13,-14)]
dG1t5 <- as.data.frame(apply(dG1t4, 2, function(x){x-mean(x)}))
names(dG1t5);colMeans(dG1t5)
# 建立二次項和交互作用項
tmp <- as.data.frame(dG1t5^2); names(tmp) <- paste0(names(tmp), "2")
dG1t <- cbind(as.data.frame(model.matrix(~(.)^2, dG1t5)), tmp)  [,-1]
names(dG1t) <- gsub(":", "", names(dG1t))
names(dG1t);colMeans(dG1t);dim(dG1t) #  spots remained
# 減去沒用的交互作用 -GbGc -GbGf -GcGd -GcGf -GcGi -GcGj -GfGi -GfGj
dG1t <- dG1t[!colnames(dG1t) %in% c("GaGc","GaGf","GbGc","GbGf","GcGd","GcGe","GcGf","GcGi","GcGj","GdGf","GeGf","GfGi","GfGj","Gc2","Gf2")]

# 生態資料
dB0t <- decostand(dt[, c(22,23,24,25)], "range", 2)
  # dB0t <- decostand(dt[, c(22,23,24,25)], 2, "standardize")
names(dB0t); dim(dB0t) # 783 spots remained

quartz(width=10, height=5)
par(cex=8/12, mar=c(5,5,0,0)+0.1, family = "Noto Sans CJK TC", mgp=c(4,1,0))
boxplot(cbind(dB0t, dG1t), xlab="變數", ylab="標準化後數值", las=2)
  # axis(1, bp, names(dG1t))
quartz.save("var-summary.pdf", type="pdf", family = "Noto Sans CJK TC")
quartz.save("var-summary.png", dpi=300, type="png", family = "Noto Sans CJK TC")
dev.off()





## 求原始生態資料降維後分數 Bscore.o
Bscore.o <- data.frame(NO = dt$NO)
fitPCA <- function(dt){
  fit <- princomp(dt, cor=F)
  score <- fit$scores[,1]
  loading <- t(loadings(fit)[,1])
  if(mean(loading < 0)) {score <- -score; loading <- -loading}
  print(summary(fit))
  print(round(loading, digits=6))
  return((score - mean(score))/sd(score))
}
fitFA <- function(dt){
  fit <- factanal(dt, 1, rotation="varimax", scores="regression", trace=T)
  score <- as.vector(fit$scores)
  loading <- t(loadings(fit))
  if(mean(loading < 0)) {score <- -score; loading <- -loading}
  fit
  print(t(loadings(fit)), digits=6)
  return((score - mean(score))/sd(score))
}
fitMDS <- function(dt){
  fit <- vegan::metaMDS(dt, distance="euclidean", k = 1, trymax = 8)
  score <- scores(fit)
  loading <- t(fit$species)
  if(mean(loading < 0)) {score <- -score; loading <- -loading}
  print(round(fit$stress, digits=6))
  print(round(loading, digits=6))
  return((score - mean(score))/sd(score))
}
fitRDA <- function(dtB, dtG){
  fit <- vegan::rda(dtB ~ . , data=dtG, scale=F)
  score <- summary(fit)$sites[,1]
  loading <- t(summary(fit)$species[,1])
  if(mean(loading < 0)) {score <- -score; loading <- -loading}
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
quartz.save("原始資料降維分數.png", dpi=300, type="png", family = "Noto Sans CJK TC")
quartz.save("原始資料降維分數.pdf", type="pdf", family = "Noto Sans CJK TC")
dev.off()







## multiple regression: diversity ~ area^1 + area^2 + interaction(among area^1)
## AIC
fit.o <- list()
fit.o$PCA.AIC <- stats::step(lm(Bscore.o$PCA ~ . , data=dG1t))
fit.o$FA.AIC  <- stats::step(lm(Bscore.o$FA  ~ . , data=dG1t))
fit.o$MDS.AIC <- stats::step(lm(Bscore.o$MDS ~ . , data=dG1t))
fit.o$RDA.AIC <- stats::step(lm(Bscore.o$RDA ~ . , data=dG1t))
fit.o$PCA.BIC <- stats::step(lm(Bscore.o$PCA ~ . , data=dG1t), k = log(nrow(dG1t)))
fit.o$FA.BIC  <- stats::step(lm(Bscore.o$FA  ~ . , data=dG1t), k = log(nrow(dG1t)))
fit.o$MDS.BIC <- stats::step(lm(Bscore.o$MDS ~ . , data=dG1t), k = log(nrow(dG1t)))
fit.o$RDA.BIC <- stats::step(lm(Bscore.o$RDA ~ . , data=dG1t), k = log(nrow(dG1t)))
# 結合係數成data.frame: coefM.o
coefM.o <- multimerge(lapply(fit.o, function(x){data.frame(coef(x))}))
  coefM.o <- merge(coefM.o, as.data.frame(coef(lm(dB0t$HsGrass ~ . , data=dG1t))), by="row.names", sort=T, all=T)
  rownames(coefM.o) <- coefM.o$Row.names
  coefM.o$`Row.names` <- NULL
  coefM.o[, ncol(coefM.o)] <- NULL
round(coefM.o, digit=3)
library(xtable)
xtable(coefM.o)

# 結合標準化係數成data.frame: coefM.o.std 使用 QuantPsyc::lm.beta()
library(QuantPsyc)
coefM.o.std <- multimerge(lapply(fit.o, function(x){data.frame(lm.beta(x))}))
  coefM.o.std <- merge(coefM.o.std, as.data.frame(coef(lm(dB0t$HsGrass ~ . , data=dG1t))), by="row.names", sort=T, all=T)
  rownames(coefM.o.std) <- coefM.o.std$Row.names
  coefM.o.std$`Row.names` <- NULL
  coefM.o.std[, ncol(coefM.o.std)] <- NULL
round(coefM.o.std, digit=3)
library(xtable)
xtable(coefM.o.std)
# % latex table generated in R 3.2.2 by xtable 1.7-4 package
# % Wed Nov  4 02:45:47 2015
# \begin{table}[ht]
# \centering
# \begin{tabular}{rrrrrrrrr}
#   \hline
#  & PCA.AIC & FA.AIC & MDS.AIC & RDA.AIC & PCA.BIC & FA.BIC & MDS.BIC & RDA.BIC \\
#   \hline
# (Intercept) &  &  &  &  &  &  &  &  \\
#   Ga & 0.15 & 0.17 & 0.07 & 0.15 & 0.16 &  &  & 0.14 \\
#   Ga2 & -0.14 & -0.20 & -0.12 &  & -0.13 &  &  &  \\
#   GaGb &  &  & 0.07 & 0.16 &  &  &  & 0.13 \\
#   GaGd & -0.34 & -0.39 & -0.35 & -0.24 & -0.30 & -0.33 & -0.31 & -0.25 \\
#   GaGe &  & -0.09 &  & 0.20 &  &  &  & 0.20 \\
#   GaGi & 0.05 &  & 0.06 & 0.10 &  &  &  & 0.09 \\
#   GaGj & -0.06 & -0.06 & -0.10 &  &  &  & -0.08 &  \\
#   Gb &  & -0.15 &  & -0.13 &  &  &  &  \\
#   Gb2 & 0.12 & 0.24 & 0.33 & 0.49 &  &  &  & 0.27 \\
#   GbGd &  &  &  &  &  &  &  &  \\
#   GbGe &  &  & 0.17 & 0.27 &  &  &  & 0.23 \\
#   GbGi & 0.09 &  & 0.11 & 0.14 &  &  &  & 0.13 \\
#   GbGj & 0.11 & 0.10 & 0.14 & 0.11 &  &  &  &  \\
#   Gc & 0.07 & 0.06 & 0.06 &  & 0.07 &  & 0.07 &  \\
#   Gd & 0.92 & 0.89 & 0.87 & 0.89 & 0.73 & 0.96 & 0.77 & 0.91 \\
#   Gd2 & -0.43 & -0.39 & -0.50 & -0.37 & -0.55 & -0.39 & -0.55 & -0.38 \\
#   GdGe & 0.32 & 0.39 & 0.24 & 0.46 &  & 0.53 &  & 0.45 \\
#   GdGi & 0.07 & 0.07 &  &  &  &  &  &  \\
#   GdGj &  &  &  &  &  &  &  &  \\
#   Ge &  &  &  &  &  &  &  &  \\
#   Ge2 & 0.07 &  & 0.11 & 0.27 &  &  &  & 0.29 \\
#   GeGi & -0.14 & -0.26 & -0.12 &  &  & -0.35 & -0.08 &  \\
#   GeGj & 0.46 & 0.42 & 0.47 & 0.45 & 0.56 & 0.43 & 0.56 & 0.46 \\
#   Gf &  &  &  &  &  &  &  &  \\
#   Gi &  &  &  &  &  &  &  &  \\
#   Gi2 &  & -0.18 &  &  &  & -0.35 &  &  \\
#   GiGj &  &  &  &  &  &  &  &  \\
#   Gj & 0.74 & 0.67 & 0.75 & 0.64 & 0.79 & 0.60 & 0.77 & 0.57 \\
#   Gj2 &  &  &  &  &  &  &  &  \\
#    \hline
# \end{tabular}
# \end{table}


## beta plot before standardizing
tmp <- data.frame(
  aiccol = 1:4,
  biccol = 5:8,
  filename = paste0("beta-", sub(".AIC", "", colnames(coefM.o)[1:4]))
  )
for(i in 1:nrow(tmp)) {
  quartz(width=7, height=3)
  par(mar=c(5,5,0,0)+0.1, family = "Noto Sans CJK TC", mgp=c(4,1,0), cex=6/12)
  plot(NULL, xlim=c(1, nrow(coefM.o)), ylim=c(-9,9), xaxt="n", yaxt="n", xlab="類型", ylab="原始迴歸係數" )
    axis(1, 1:nrow(coefM.o), rownames(coefM.o), las=2)
    axis(2, seq(-10, 10, 2), c( -(2^(abs(seq(-10, -2, 2)))), 0, 2^(seq(2, 10, 2)) ), las=2)
    abline(h=0)
    arrows(
      1:nrow(coefM.o)-0.2,
      rep(0, nrow(coefM.o)),
      1:nrow(coefM.o)-0.2,
      ifelse(
        coefM.o[, tmp$aiccol[i]]>=0|is.na(coefM.o[, tmp$aiccol[i]]),
        log(coefM.o[,tmp$aiccol[i]]+1, 2),
        -log(-coefM.o[,tmp$aiccol[i]]+1, 2)
      ),
      code=2, length = 0.05, col=1 #, lwd=abs(coefM.o[,1])/max(abs(coefM.o[,1]))*5+1
    )
    arrows(
      1:nrow(coefM.o)+0.2,
      rep(0, nrow(coefM.o)),
      1:nrow(coefM.o)+0.2,
      ifelse(
        coefM.o[,tmp$biccol[i]]>=0|is.na(coefM.o[,tmp$biccol[i]]),
        log(coefM.o[,tmp$biccol[i]]+1, 2),
        -log(-coefM.o[,tmp$biccol[i]]+1, 2)
      ) ,
      code=2, length = 0.05, col=2, lty=1#, lwd=abs(coefM.o[,2])/max(abs(coefM.o[,2]))*5+1
    )
    legend(1, 8, c("AIC", "BIC"), lty=c(1,1), col=c(1,2))
  quartz.save( paste0(tmp$filename[i], ".png") , dpi=300, type="png", family = "Noto Sans CJK TC")
  quartz.save( paste0(tmp$filename[i], ".pdf"), type="pdf", family = "Noto Sans CJK TC")
  dev.off()
}


## beta plot after standardizing
tmp <- data.frame(
  aiccol = 1:4,
  biccol = 5:8,
  filename = paste0("beta-std", sub(".AIC", "", colnames(coefM.o.std)[1:4]))
)
for(i in 1:nrow(tmp)) {
  quartz(width=7, height=3)
  par(mar=c(5,5,0,0)+0.1, family = "Noto Sans CJK TC", mgp=c(4,1,0), cex=6/12)
  plot(NULL, xlim=c(1, nrow(coefM.o.std)), ylim=c(-3.5,3.5), xaxt="n", yaxt="n", xlab="類型", ylab="標準化迴歸係數" )
    axis(1, 1:nrow(coefM.o.std), rownames(coefM.o.std), las=2)
    axis(2, seq(-3.5, 3.5, 0.5), seq(-3.5, 3.5, 0.5), las=2)
    abline(h=0)
    arrows(
      1:nrow(coefM.o.std)-0.2,
      rep(0, nrow(coefM.o.std)),
      1:nrow(coefM.o.std)-0.2,
      coefM.o.std[,tmp$aiccol[i]],
      code=2, length = 0.05, col=1 #, lwd=abs(coefM.o[,1])/max(abs(coefM.o[,1]))*5+1
    )
    arrows(
      1:nrow(coefM.o)+0.2,
      rep(0, nrow(coefM.o)),
      1:nrow(coefM.o)+0.2,
      coefM.o.std[,tmp$biccol[i]],
      code=2, length = 0.05, col=2, lty=1#, lwd=abs(coefM.o[,2])/max(abs(coefM.o[,2]))*5+1
    )
    legend(1, 2, c("AIC", "BIC"), lty=c(1,1), col=c(1,2))
  quartz.save( paste0(tmp$filename[i], ".png") , dpi=300, type="png", family = "Noto Sans CJK TC")
  quartz.save( paste0(tmp$filename[i], ".pdf"), type="pdf", family = "Noto Sans CJK TC")
  dev.off()
}

## beta plot after standardizing


########### find best BAF
dB0t.pred <- data.frame(
  # `CCA.AIC.p` = predict(lm.cca.aic, dG1t),
  # `CCA.BIC.p` = predict(lm.cca.bic, dG1t),
  MDS.AIC.p = predict(fit.o$MDS.AIC, dG1t),
  MDS.BIC.p = predict(fit.o$MDS.BIC, dG1t),
  FA.AIC.p  = predict(fit.o$FA.AIC, dG1t),
  FA.BIC.p  = predict(fit.o$FA.BIC, dG1t),
  PCA.AIC.p = predict(fit.o$PCA.AIC, dG1t),
  PCA.BIC.p = predict(fit.o$PCA.BIC, dG1t),
  RDA.AIC.p = predict(fit.o$RDA.AIC, dG1t),
  RDA.BIC.p = predict(fit.o$RDA.BIC, dG1t)
  # RDA.AICd.p = predict(mlm.rda.bic, dG1t),
  # RDA.BICd.p = predict(mlm.rda.bic, dG1t)
)

quartz(width=9, height=9)
par(mfrow=c(1,4), cex=6/12)
pairs(
  cbind(dB0t.pred, dB0t),
  lower.panel=panel.smooth,
  upper.panel=panel.cor, pch=".")
quartz.save("all.png", dpi=300, type="png", family = "Noto Sans CJK TC")
quartz.save("all.pdf", type="pdf", family = "Noto Sans CJK TC")
dev.off()















# 創造 10-fold 並回傳 trainning set
# 切出 TD & VD
library(caret)
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
Bscore.VD <- list()
fit.TD <- list()
Bscore.TD.pred <- list()
for(i in 1:length(TD.row)){
# for(i in 1:2){
  # 結合已降維的 dB0t.TD[[i]] 成 Bscore.TD[[i]]
  Bscore.TD[[i]] <- data.frame(NO = dt$NO[TD.row[[i]]])
  Bscore.TD[[i]]$PCA <- Bscore.o$PCA[TD.row[[i]]]
  Bscore.TD[[i]]$FA  <- Bscore.o$FA[TD.row[[i]]]
  Bscore.TD[[i]]$MDS <- Bscore.o$MDS[TD.row[[i]]]
  Bscore.TD[[i]]$RDA <- Bscore.o$RDA[TD.row[[i]]]
  rownames(Bscore.TD[[i]]) <- dt$NO[TD.row[[i]]]
  Bscore.TD[[i]]$NO <- NULL

  # 結合已降維的 dB0t.VD[[i]] 成 Bscore.VD[[i]]
  Bscore.VD[[i]] <- data.frame(NO = dt$NO[VD.row[[i]]])
  Bscore.VD[[i]]$PCA.AIC <- Bscore.VD[[i]]$PCA.BIC <- Bscore.o$PCA[VD.row[[i]]]
  Bscore.VD[[i]]$FA.AIC  <- Bscore.VD[[i]]$FA.BIC  <- Bscore.o$FA[VD.row[[i]]]
  Bscore.VD[[i]]$MDS.AIC <- Bscore.VD[[i]]$MDS.BIC <- Bscore.o$MDS[VD.row[[i]]]
  Bscore.VD[[i]]$RDA.AIC <- Bscore.VD[[i]]$RDA.BIC <- Bscore.o$RDA[VD.row[[i]]]
  rownames(Bscore.VD[[i]]) <- dt$NO[VD.row[[i]]]
  Bscore.VD[[i]]$NO <- NULL

  # multiple regression: Bscore.TD[[i]]$xxx ~ ., data=dG1t.TD
  fit.TD[[i]] <- list()
  fit.TD[[i]]$PCA.AIC <- stats::step(lm(Bscore.TD[[i]]$PCA ~ . , data=dG1t.TD[[i]]))
  fit.TD[[i]]$PCA.BIC <- stats::step(lm(Bscore.TD[[i]]$PCA ~ . , data=dG1t.TD[[i]]), k = log(nrow(dG1t.TD[[i]])))
  fit.TD[[i]]$FA.AIC  <- stats::step(lm(Bscore.TD[[i]]$FA  ~ . , data=dG1t.TD[[i]]))
  fit.TD[[i]]$FA.BIC  <- stats::step(lm(Bscore.TD[[i]]$FA  ~ . , data=dG1t.TD[[i]]), k = log(nrow(dG1t.TD[[i]])))
  fit.TD[[i]]$MDS.AIC <- stats::step(lm(Bscore.TD[[i]]$MDS ~ . , data=dG1t.TD[[i]]))
  fit.TD[[i]]$MDS.BIC <- stats::step(lm(Bscore.TD[[i]]$MDS ~ . , data=dG1t.TD[[i]]), k = log(nrow(dG1t.TD[[i]])))
  fit.TD[[i]]$RDA.AIC <- stats::step(lm(Bscore.TD[[i]]$RDA ~ . , data=dG1t.TD[[i]]))
  fit.TD[[i]]$RDA.BIC <- stats::step(lm(Bscore.TD[[i]]$RDA ~ . , data=dG1t.TD[[i]]), k = log(nrow(dG1t.TD[[i]])))

  # 以訓練出的 fit.TD$xxx.yIC 模型，求 dG1t.VD 的預測值
  Bscore.TD.pred[[i]] <- data.frame(
    PCA.AIC.p = predict(fit.TD[[i]]$PCA.AIC, dG1t.VD[[i]]),
    PCA.BIC.p = predict(fit.TD[[i]]$PCA.BIC, dG1t.VD[[i]]),
    FA.AIC.p  = predict(fit.TD[[i]]$FA.AIC, dG1t.VD[[i]]),
    FA.BIC.p  = predict(fit.TD[[i]]$FA.BIC, dG1t.VD[[i]]),
    MDS.AIC.p = predict(fit.TD[[i]]$MDS.AIC, dG1t.VD[[i]]),
    MDS.BIC.p = predict(fit.TD[[i]]$MDS.BIC, dG1t.VD[[i]]),
    RDA.AIC.p = predict(fit.TD[[i]]$RDA.AIC, dG1t.VD[[i]]),
    RDA.BIC.p = predict(fit.TD[[i]]$RDA.BIC, dG1t.VD[[i]])
  )
}
# 結合所有預測 Bscore.TD.pred.final 和所有驗證 Bscore.VD.final
Bscore.TD.pred.求最佳迴歸法  <- do.call("rbind", Bscore.TD.pred)
  Bscore.TD.pred.求最佳迴歸法 <- Bscore.TD.pred.求最佳迴歸法[order(as.numeric(rownames(Bscore.TD.pred.求最佳迴歸法))), ]
Bscore.VD.求最佳迴歸法 <- do.call("rbind", Bscore.VD)
  Bscore.VD.求最佳迴歸法 <- Bscore.VD.求最佳迴歸法[order(as.numeric(rownames(Bscore.VD.求最佳迴歸法))), ]
dim(Bscore.TD.pred.求最佳迴歸法) ; dim(Bscore.VD.求最佳迴歸法)
names(Bscore.TD.pred.求最佳迴歸法) ; names(Bscore.VD.求最佳迴歸法)
# 看看誰最好
library(hydroGOF)
(tmp <- rbind(
  RMSE = rmse(Bscore.TD.pred.求最佳迴歸法, Bscore.VD.求最佳迴歸法),
  NRMSE = nrmse(Bscore.TD.pred.求最佳迴歸法, Bscore.VD.求最佳迴歸法)/100
))
library(xtable)
xtable(tmp)
# seed = 52004800
# PCA.AIC.p PCA.BIC.p FA.AIC.p  FA.BIC.p MDS.AIC.p MDS.BIC.p RDA.AIC.p RDA.BIC.p
# RMSE   3.084222  3.082647 2.648136 0.8098411  3.703313  3.947621  2.853432   2.95322
# NRMSE  3.084000  3.083000 2.648000 0.8100000  3.703000  3.948000  2.853000   2.95300
# \begin{table}[ht]
# \centering
# \begin{tabular}{rrrrrrrrr}
#   \hline
#  & PCA.AIC.p & PCA.BIC.p & FA.AIC.p & FA.BIC.p & MDS.AIC.p & MDS.BIC.p & RDA.AIC.p & RDA.BIC.p \\
#   \hline
# RMSE & 3.08 & 3.08 & 2.65 & 0.81 & 3.70 & 3.95 & 2.85 & 2.95 \\
#   NRMSE & 3.08 & 3.08 & 2.65 & 0.81 & 3.70 & 3.95 & 2.85 & 2.95 \\
#    \hline
# \end{tabular}
# \end{table}















#### 10-fold CV 找8種求解法何者較佳
#### 目標：求最佳整套求解方式（包括降維與迴歸方式）
Bscore.TD <- list()
Bscore.VD <- list()
fit.TD <- list()
Bscore.TD.pred <- list()
for(i in 1:length(TD.row)){
  # 生態訓練資料 dB0t.TD[[i]] 重新降維成 Bscore.TD[[i]]
  Bscore.TD[[i]] <- data.frame(NO = dt$NO[TD.row[[i]]])
  Bscore.TD[[i]]$PCA <- fitPCA(dB0t.TD[[i]])
  Bscore.TD[[i]]$FA  <- fitFA(dB0t.TD[[i]])
  Bscore.TD[[i]]$MDS <- fitMDS(dB0t.TD[[i]])
  Bscore.TD[[i]]$RDA <- fitRDA(dB0t.TD[[i]], dG1t.TD[[i]])
  rownames(Bscore.TD[[i]]) <- dt$NO[TD.row[[i]]]
  Bscore.TD[[i]]$NO <- NULL

  # 生態驗證資料 dB0t.VD[[i]] = Bscore.VD[[i]]
  Bscore.VD[[i]] <- data.frame(NO = dt$NO[VD.row[[i]]])
  Bscore.VD[[i]]$PCA.AIC <- Bscore.VD[[i]]$PCA.BIC <- fitPCA(dB0t.VD[[i]])
  Bscore.VD[[i]]$FA.AIC  <- Bscore.VD[[i]]$FA.BIC  <- fitFA(dB0t.VD[[i]])
  Bscore.VD[[i]]$MDS.AIC <- Bscore.VD[[i]]$MDS.BIC <- fitMDS(dB0t.VD[[i]])
  Bscore.VD[[i]]$RDA.AIC <- Bscore.VD[[i]]$RDA.BIC <- fitRDA(dB0t.VD[[i]], dG1t.VD[[i]])
  rownames(Bscore.VD[[i]]) <- dt$NO[VD.row[[i]]]
  Bscore.VD[[i]]$NO <- NULL

  # multiple regression: Bscore.TD[[i]]$xxx ~ ., data=dG1t.TD
  fit.TD[[i]] <- list()
  fit.TD[[i]]$PCA.AIC <- stats::step(lm(Bscore.TD[[i]]$PCA ~ . , data=dG1t.TD[[i]]))
  fit.TD[[i]]$PCA.BIC <- stats::step(lm(Bscore.TD[[i]]$PCA ~ . , data=dG1t.TD[[i]]), k = log(nrow(dG1t.TD[[i]])))
  fit.TD[[i]]$FA.AIC  <- stats::step(lm(Bscore.TD[[i]]$FA  ~ . , data=dG1t.TD[[i]]))
  fit.TD[[i]]$FA.BIC  <- stats::step(lm(Bscore.TD[[i]]$FA  ~ . , data=dG1t.TD[[i]]), k = log(nrow(dG1t.TD[[i]])))
  fit.TD[[i]]$MDS.AIC <- stats::step(lm(Bscore.TD[[i]]$MDS ~ . , data=dG1t.TD[[i]]))
  fit.TD[[i]]$MDS.BIC <- stats::step(lm(Bscore.TD[[i]]$MDS ~ . , data=dG1t.TD[[i]]), k = log(nrow(dG1t.TD[[i]])))
  fit.TD[[i]]$RDA.AIC <- stats::step(lm(Bscore.TD[[i]]$RDA ~ . , data=dG1t.TD[[i]]))
  fit.TD[[i]]$RDA.BIC <- stats::step(lm(Bscore.TD[[i]]$RDA ~ . , data=dG1t.TD[[i]]), k = log(nrow(dG1t.TD[[i]])))

  # 以 fit.TD[[i]]$xxx.yIC 模型，求 dG1t.TD 的預測值 Bscore.TD.pred
  Bscore.TD.pred[[i]] <- data.frame(
    PCA.AIC.p = predict(fit.TD[[i]]$PCA.AIC, dG1t.VD[[i]]),
    PCA.BIC.p = predict(fit.TD[[i]]$PCA.BIC, dG1t.VD[[i]]),
    FA.AIC.p  = predict(fit.TD[[i]]$FA.AIC, dG1t.VD[[i]]),
    FA.BIC.p  = predict(fit.TD[[i]]$FA.BIC, dG1t.VD[[i]]),
    MDS.AIC.p = predict(fit.TD[[i]]$MDS.AIC, dG1t.VD[[i]]),
    MDS.BIC.p = predict(fit.TD[[i]]$MDS.BIC, dG1t.VD[[i]]),
    RDA.AIC.p = predict(fit.TD[[i]]$RDA.AIC, dG1t.VD[[i]]),
    RDA.BIC.p = predict(fit.TD[[i]]$RDA.BIC, dG1t.VD[[i]])
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
# 看看誰最好
library(hydroGOF)
(tmp <- rbind(
  RMSE = rmse(Bscore.TD.pred.final, Bscore.VD.final),
  NRMSE = nrmse(Bscore.TD.pred.final, Bscore.VD.final)/100
))
(tmp <- rbind(
  RMSE = rmse(Bscore.TD.pred.final[-418,], Bscore.VD.final[-418,]),
  NRMSE = nrmse(Bscore.TD.pred.final[-418,], Bscore.VD.final[-418,])/100
))
library(xtable)
xtable(tmp)
# seed = 52004800
# PCA.AIC.p PCA.BIC.p FA.AIC.p FA.BIC.p    NMDS1    NMDS1 RDA.AIC.p RDA.BIC.p
# RMSE   3.045977   3.04736 2.638458 0.875508 3.670323 3.932308  2.770957   2.88824
# NRMSE  3.063000   3.06400 2.653000 0.880000 3.691000 3.954000  2.786000   2.90400
# \begin{table}[ht]
# \centering
# \begin{tabular}{rrrrrrrrr}
#   \hline
#  & PCA.AIC.p & PCA.BIC.p & FA.AIC.p & FA.BIC.p & NMDS1 & NMDS1 & RDA.AIC.p & RDA.BIC.p \\
#   \hline
# RMSE & 3.05 & 3.05 & 2.64 & 0.88 & 3.67 & 3.93 & 2.77 & 2.89 \\
#   NRMSE & 3.06 & 3.06 & 2.65 & 0.88 & 3.69 & 3.95 & 2.79 & 2.90 \\
#    \hline
# \end{tabular}
# \end{table}

# seed = 19821006
# PCA.AIC.p PCA.BIC.p FA.AIC.p  FA.BIC.p    NMDS1    NMDS1 RDA.AIC.p RDA.BIC.p
# RMSE   1.054052  1.054443 2.628768 0.8632065 1.086618 1.162657  0.521941 0.5433606
# NRMSE  3.060000  3.061000 2.689000 0.8830000 3.158000 3.379000  1.949000 2.0290000

# mean squared error (MSE)
# root-mean-square deviation (RMSD)
  # sqrt(colSums((Bscore.TD.pred.final - Bscore.VD.final)^2 / nrow(Bscore.TD.pred.final)))
# median absolute deviation (MAD)
# mse(Bscore.TD.pred.final, Bscore.VD.final)
# library(reshape2)
# tmp <- (Bscore.TD.pred.final - Bscore.VD.final)^2
# summary(tmp)
boxplot(
  apply( sqrt((Bscore.TD.pred.final - Bscore.VD.final)^2/nrow(Bscore.TD.pred.final)) , 2 , function(x){return(x/sd(x))} )
)
boxplot(
  apply( sqrt((Bscore.TD.pred.final[-418, ] - Bscore.VD.final[-418, ])^2/nrow(Bscore.TD.pred.final[-418, ])) , 2 , function(x){return(x/sd(x))} )
)
max(Bscore.TD.pred.final[,1])
dt [max(Bscore.TD.pred.final[,1]) == Bscore.TD.pred.final[,1] , ]
# kruskal.test(value ~ Var2,  data=melt(tmp))
# library(PMCMR)
# posthoc.kruskal.nemenyi.test(melt(tmp)$value, melt(tmp)$Var2, method="Tukey")















#### 10-fold CV 找哪種降維方法最好 ---- 不可能，因為無法訓練出 FA MDS 的模型來預測分數
#### 改做 bootstrap kendall correlation 找哪種降維方法最好
k <- 10000
ind <- vector("list", k)
set.seed(52004800)#; set.seed(19821006); set.seed(124234545)
ind <- lapply(ind, function(x){return(sample(1:nrow(dB0t), replace=T))})
bootR <- vector("list", k)
pb <- txtProgressBar(min = 0, max = k, style = 3)
for(i in 1:k){
  bootR[[i]] <- data.frame(
    cor(Bscore.o[ind[[i]],-1]$PCA, dB0t[ind[[i]],], method="kendall" ),
    cor(Bscore.o[ind[[i]],-1]$FA,  dB0t[ind[[i]],], method="kendall" ),
    cor(Bscore.o[ind[[i]],-1]$MDS, dB0t[ind[[i]],], method="kendall" ),
    cor(Bscore.o[ind[[i]],-1]$RDA, dB0t[ind[[i]],], method="kendall" )
  )
  setTxtProgressBar(pb, i)
}
close(pb)
bootR.final <- do.call("rbind", bootR)
colnames(bootR.final) <- paste0(
  rep(c("HsWood", "HsGrass", "HsSpider", "HsInsect"), 4),
  c(rep(".PCA", 4), rep(".FA", 4), rep(".MDS", 4), rep(".RDA", 4))
)
apply(bootR.final, 2, function(x){
    quantile(x, c(0.0005,0.005,0.025,0.975,0.995,0.9995))
})
save(bootR.final, file="boorT.final.Rdata")
