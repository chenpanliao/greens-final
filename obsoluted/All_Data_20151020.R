panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    ra <- abs(cor(x, y))
    r <- cor(x, y)
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
BAF.mt <- matrix(
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
# 去除 Gx > 9000
# 面積以列標準化為 sum = 1
library(vegan)
dt <- na.omit(d); dim(dt)
dt <- subset(dt, Gx+Ga+Gb+Gc+Gd+Ge+Gf+Gg+Gh+Gi+Gj+Gk+Gl+Gm >= 7000) ; dim(dt)
dt <- subset(dt, Gx+Ga+Gb+Gc+Gd+Ge+Gf+Gg+Gh+Gi+Gj+Gk+Gl+Gm <= 13000) ; dim(dt)
dt <- subset(dt, Gx <= 9000) ; dim(dt)
dt <- subset(dt, Gx+Ga+Gb+Gc+Gd+Ge+Gf+Gg+Gh+Gi+Gj+Gk+Gl+Gm >= 7000 & Gx+Ga+Gb+Gc+Gd+Ge+Gf+Gg+Gh+Gi+Gj+Gk+Gl+Gm <= 13000 & Gx <= 9000)
dG0t1 <- dt[, c(4:17)]
dG0t2 <- dG0t1[, c(-8, -12,-13,-14)]
names(dG0t2)
dG0t3 <- decostand(dG0t2, "total", 1)
names(dG0t3);rowSums(dG0t3)
dG1t4 <- dG0t3[, c(-1, -8, -12,-13,-14)]
dG1t5 <- as.data.frame(apply(dG1t4, 2, function(x){x-mean(x)}))
names(dG1t5);colMeans(dG1t5)

tmp <- as.data.frame(dG1t5^2); names(tmp) <- paste0(names(tmp), "2")
dG1t <- cbind(as.data.frame(model.matrix(~(.)^2, dG1t5)), tmp)  [,-1]
names(dG1t) <- gsub(":", "", names(dG1t))
names(dG1t);colMeans(dG1t);dim(dG1t) #  spots remained

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
# library(tidyr)
# data_long <- gather(cbind(dB0t,dG1t), var, val)
# data_long
# names(data_long)
# ggplot(data_long, aes(factor(var), val)) + geom_violin()




## 求生物 MDS 分數 Bscore.mds
mds1 <- vegan::metaMDS(dB0t, distance="euclidean", k=1, trymax = 50)
Bscore.mds <- -scores(mds1)
round(t(-mds1$species), digits=6)
mds1
# Stress:     0.2842971
#       HsWood   HsGrass HsSpider HsInsect
# MDS1 0.73107 -0.009557 0.683752 0.207176
tm <- cbind(`MDS score`=Bscore.mds[,1], dB0t)
quartz(width=7, height=2)
par(mfrow=c(1,4), cex=6/12)
for(i in 2:ncol(tm)){
  plot(tm[[1]] ~ tm[[i]], xlab=names(tm)[i], ylab="MDS1 socre", pch=".")
  lines(lowess(tm[[i]], tm[[1]], f=4/5, iter = 100), col = 2)
  title(paste0("r = ", round(cor(tm[[1]], tm[[i]]), digits=2)))
}
quartz.save("scoreMDS.pdf", type="pdf", family = "Noto Sans CJK TC")
quartz.save("scoreMDS.png", dpi=300, type="png", family = "Noto Sans CJK TC")
dev.off()


## 求生物 PCA 分數 Bsocre.pca
pca1 <- princomp(dB0t, cor=F)
Bscore.pca <- -pca1$scores
summary(pca1)
round(t(-loadings(pca1)[,1]), digits=6)
tm <- cbind(Bscore.pca[,1], dB0t)
quartz(width=7, height=2)
par(mfrow=c(1,4), cex=6/12)
for(i in 2:ncol(tm)){
  plot(tm[[1]] ~ tm[[i]], xlab=names(tm)[i], ylab="PC1 socre", pch=".")
  lines(lowess(tm[[i]], tm[[1]], f=4/5, iter = 100), col = 2)
  title(paste0("r = ", round(cor(tm[[1]], tm[[i]]), digits=2)))
}
quartz.save("scorePCA.pdf", type="pdf", family = "Noto Sans CJK TC")
quartz.save("scorePCA.png", dpi=300, type="png", family = "Noto Sans CJK TC")
dev.off()


## 求 RDA1 分數對應的係數
library(vegan)
g <- paste0(c(
  "cbind(", paste0(names(dB0t), sep="", collapse=","), ")"
  ), collapse="")
h <- paste0(
  c(
    paste0(names(dG1t), "", collapse="+")
  ), collapse=" + "
)
as.formula(paste(g, "~", h))
rda1 <- vegan::rda(
  # as.formula(paste("dB0t ~", h)),
  dB0t ~ . ,
  data=dG1t, scale=F)
round(t(-summary(rda1)$species[,1]), digits=6)
head(summary(rda1))
Bscore.rda <- -summary(rda1)$sites
tm <- cbind(Bscore.rda[,1], dB0t)
quartz(width=7, height=2)
par(mfrow=c(1,4), cex=6/12)
for(i in 2:ncol(tm)){
  plot(tm[[1]] ~ tm[[i]], xlab=names(tm)[i], ylab="RD1 socre", pch=".")
  lines(lowess(tm[[i]], tm[[1]], f=4/5, iter = 100), col = 2)
  title(paste0("r = ", round(cor(tm[[1]], tm[[i]]), digits=2)))
}
quartz.save("scoreRDA.png", dpi=300, type="png", family = "Noto Sans CJK TC")
quartz.save("scoreRDA.pdf", type="pdf", family = "Noto Sans CJK TC")
dev.off()




## 求生物 CCA 分數 Bsocre.cca
g <- paste0(c(
  "cbind(", paste0(names(dB0t), sep="", collapse=","), ")"
  ), collapse="")
h <- paste0(
  c(
    paste0(names(dG1t), "", collapse="+"),
    paste0("I(", names(dG1t), "^2)", sep="", collapse="+"),
    paste0(c("(",paste0(names(dG1t), "", collapse="+"),")^2"), collapse="")
  ), collapse=" + "
)
as.formula(paste(g, "~", h))
cc1 <- candisc::cancor(cbind(HsWood, HsGrass, HsSpider, HsInsect) ~ Ga + Gb + Gc + Gd +
    Ge + Gf + Gi + Gj + I(Ga^2) + I(Gb^2) + I(Gc^2) + I(Gd^2) +
    I(Ge^2) + I(Gf^2) + I(Gi^2) + I(Gj^2) + (Ga + Gb + Gc + Gd +
    Ge + Gf + Gi + Gj)^2 -Gc:Gf,
    data=cbind(dB0t,dG1t))
    ## 預先拿到很多 interaction 不然太多 0
# cc1 <- candisc::cancor(cbind(HsWood, HsGrass, HsSpider, HsInsect) ~ (Ga + Gb + Gc + Gd +
#     Ge + Gf + Gi + Gj)^2,
#     data=cbind(dB0t,dG1t))
summary(cc1)
Bscore.cca <- cc1$scores$X
round(t(cc1$coef$Y[,1]), digits=6)
tm <- cbind(`CCA score`=Bscore.cca[,1], dB0t)
quartz(width=7, height=2)
par(mfrow=c(1,4), cex=6/12)
for(i in 2:ncol(tm)){
  plot(tm[[1]] ~ tm[[i]], xlab=names(tm)[i], ylab="CC1 socre", pch=".")
  lines(lowess(tm[[i]], tm[[1]], f=4/5, iter = 100), col = 2)
  title(paste0("r = ", round(cor(tm[[1]], tm[[i]]), digits=2)))
}
quartz.save("scoreCCA.pdf", type="pdf", family = "Noto Sans CJK TC")
quartz.save("scoreCCA.png", dpi=300, type="png", family = "Noto Sans CJK TC")
dev.off()


## 求 FA1 分數對應的係數
fa1 <- factanal(dB0t, 1, rotation="varimax", scores="regression", trace=T)
print(t(fa1$loadings), digits=6)
fa1
Bscore.fa <- data.frame(FA1 = fa1$scores)
tm <- cbind(Bscore.fa, dB0t)
quartz(width=7, height=2)
par(mfrow=c(1,4), cex=6/12)
for(i in 2:ncol(tm)){
  plot(tm[[1]] ~ tm[[i]], xlab=names(tm)[i], ylab="FA(1) socre", pch=".")
  lines(lowess(tm[[i]], tm[[1]], f=4/5, iter = 100), col = 2)
  title(paste0("r = ", round(cor(tm[[1]], tm[[i]]), digits=2)))
}
quartz.save("scoreFA.png", dpi=300, type="png", family = "Noto Sans CJK TC")
quartz.save("scoreFA.pdf", type="pdf", family = "Noto Sans CJK TC")
dev.off()











## multiple regression: diversity ~ area^1 + area^2 + interaction(among area^1)

## AIC
lm.cca.aic <- stats::step(lm(Bscore.cca[,1] ~ . , data=dG1t))
  lm.cca.aic.coef <- as.data.frame(coef(lm.cca.aic))
lm.fa.aic <- stats::step(lm(Bscore.fa[,1] ~ . , data=dG1t))
  lm.fa.aic.coef <- as.data.frame(coef(lm.fa.aic))
lm.mds.aic <- stats::step(lm(Bscore.mds[,1] ~ . , data=dG1t))
  lm.mds.aic.coef <- as.data.frame(coef(lm.mds.aic))
lm.pca.aic <- stats::step(lm(Bscore.pca[,1] ~ . , data=dG1t))
  lm.pca.aic.coef <- as.data.frame(coef(lm.pca.aic))
lm.rda.aic <- stats::step(lm(Bscore.rda[,1] ~ . , data=dG1t))
  lm.rda.aic.coef <- as.data.frame(coef(lm.rda.aic))
## BIC
lm.cca.bic <- stats::step(lm(Bscore.cca[,1] ~ . , data=dG1t), k = log(nrow(dB0t)))
  lm.cca.bic.coef <- as.data.frame(coef(lm.cca.bic))
lm.fa.bic <- stats::step(lm(Bscore.fa[,1] ~ . , data=dG1t), k = log(nrow(dB0t)))
  lm.fa.bic.coef <- as.data.frame(coef(lm.fa.bic))
lm.mds.bic <- stats::step(lm(Bscore.mds[,1] ~ . , data=dG1t), k = log(nrow(dB0t)))
  lm.mds.bic.coef <- as.data.frame(coef(lm.mds.bic))
lm.pca.bic <- stats::step(lm(Bscore.pca[,1] ~ . , data=dG1t), k = log(nrow(dB0t)))
  lm.pca.bic.coef <- as.data.frame(coef(lm.pca.bic))
lm.rda.bic <- stats::step(lm(Bscore.rda[,1] ~ . , data=dG1t), k = log(nrow(dB0t)))
  lm.rda.bic.coef <- as.data.frame(coef(lm.rda.bic))


multimerge <- function(x, ...) {
  if(!is.list(x)) stop("x must be a list.")
  y <- x[[1]]
  i <- 2
  while(i <= length(x)){
    y <- merge(y, x[[i]], by='row.names', all=T)
    rownames(y) <- y$`Row.names`
    y$`Row.names` <- NULL
    i <- i + 1
  }
  names(y) <- names(x)
  return(y)
}
coefM <- multimerge(list(
  # CCA.AIC=lm.cca.aic.coef,
  # CCA.BIC=lm.cca.bic.coef,
  MDS.AIC=lm.mds.aic.coef,
  MDS.BIC=lm.mds.bic.coef,
  PCA.AIC=lm.pca.aic.coef,
  PCA.BIC=lm.pca.bic.coef,
  RDA.AIC=lm.rda.aic.coef,
  RDA.BIC=lm.rda.bic.coef,
  FA.AIC=lm.fa.aic.coef,
  FA.BIC=lm.fa.bic.coef)) # )[-1,]
# rownames(coefM) <- gsub("`", "", rownames(coefM))
# rownames(coefM) <- gsub(":", "", rownames(coefM))
round(coefM, digits=3)[,]
# TBAF <- decostand(coefM, "range", 2, na.rm=T)*10
# round(TBAF, digits=4)

## beta plot CCA

tmp <- data.frame(
  aiccol = seq(1,7,2),
  biccol = seq(2,8,2),
  filename = paste0("beta-", sub(".AIC", "", colnames(coefM)[seq(1,7,2)])),
  ylabname = paste0("log_2 (", sub(".AIC", "", colnames(coefM)[seq(1,8,2)]), " 迴歸係數)")
  )
for(i in 1:nrow(tmp)) {
  quartz(width=7, height=3)
  par(mar=c(5,5,0,0)+0.1, family = "Noto Sans CJK TC", mgp=c(4,1,0), cex=9/12)
  plot(NULL, xlim=c(1, nrow(coefM)), ylim=c(-9,9), xaxt="n", yaxt="n", xlab="類型", ylab=tmp$ylabname[i] )
    axis(1, 1:nrow(coefM), rownames(coefM), las=2)
    axis(2, seq(-10, 10, 2), c( -(2^(abs(seq(-10, -2, 2)))), 0, 2^(seq(2, 10, 2)) ), las=2)
    abline(h=0)
    arrows(
      1:nrow(coefM)-0.2,
      rep(0, nrow(coefM)),
      1:nrow(coefM)-0.2,
      ifelse(coefM[, tmp$aiccol[i]]>=0, log(coefM[,tmp$aiccol[i]]+1, 2), -log(-coefM[,tmp$aiccol[i]]+1, 2)) ,
      code=2, length = 0.05, lend=3, col=1 #, lwd=abs(coefM[,1])/max(abs(coefM[,1]))*5+1
    )
    arrows(
      1:nrow(coefM)+0.2,
      rep(0, nrow(coefM)),
      1:nrow(coefM)+0.2,
      ifelse(coefM[,tmp$biccol[i]]>=0, log(coefM[,tmp$biccol[i]]+1, 2), -log(-coefM[,tmp$biccol[i]]+1, 2)) ,
      code=2, length = 0.05, lend=1, col=2, lty=2#, lwd=abs(coefM[,2])/max(abs(coefM[,2]))*5+1
    )
    legend(1, 8, c("AIC", "BIC"), lty=c(1,2), col=c(1,2))
  quartz.save( paste0(tmp$filename[i], ".png") , dpi=300, type="png", family = "Noto Sans CJK TC")
  quartz.save( paste0(tmp$filename[i], ".pdf"), type="pdf", family = "Noto Sans CJK TC")
  dev.off()
}




########### find best BAF
dB0t.pred <- data.frame(
  # `CCA.AIC.p` = predict(lm.cca.aic, dG1t),
  # `CCA.BIC.p` = predict(lm.cca.bic, dG1t),
  MDS.AIC.p = predict(lm.mds.aic, dG1t),
  MDS.BIC.p = predict(lm.mds.bic, dG1t),
  FA.AIC.p = predict(lm.fa.aic, dG1t),
  FA.BIC.p = predict(lm.fa.bic, dG1t),
  PCA.AIC.p = predict(lm.pca.aic, dG1t),
  PCA.BIC.p = predict(lm.pca.bic, dG1t),
  RDA.AIC.p = predict(lm.rda.aic, dG1t),
  RDA.BIC.p = predict(lm.rda.bic, dG1t)
  # RDA.AICd.p = predict(mlm.rda.bic, dG1t),
  # RDA.BICd.p = predict(mlm.rda.bic, dG1t)
)

quartz(width=7, height=7)
par(mfrow=c(1,4), cex=4/12)
pairs(
  cbind(dB0t.pred, dB0t),
  lower.panel=panel.smooth,
  upper.panel=panel.cor, pch=".")
quartz.save("all.png", dpi=300, type="png", family = "Noto Sans CJK TC")
quartz.save("all.pdf", type="pdf", family = "Noto Sans CJK TC")
dev.off()

















# tc <- trainControl(method = "repeatedcv", number = 10, repeats = 1)
# (lm.mds.aic.train <- train(Bscore.mds[,1] ~ . , data=dG1t, method='lm', trControl = tc))
# Root Mean Square Error
# lm.mds.aic.train$resample


# 創造 10-fold 並回傳 trainning set
library(caret)
createFolds(1:100, k = 10, list = TRUE, returnTrain = F)
createFolds(1:100, k = 10, list = TRUE, returnTrain = F)

#### 10-fold CV 找6種求解法何者較佳
# 選方案 6 選 1
# 1. 切出 90% td 與 10% vd
# 2. 把 td 中生態資料降維，與環境資料做 multiple regression，得到訓練出的 lm.td
# 3. 以 vd 環境資料丟入 lm.td 中，取得預期值
# 4. repeat 1-3 共 10 次，得到所有樣本的生態預期值
# 5. 所有生態實際值 - 所有生態預期值 的 mean 和 sd 為所求

####
# 切出 TD & VD
ad.row <- rep(list(1:nrow(dB0t)), 10)
set.seed(52004800)
TD.row <- createFolds(1:nrow(dB0t), k = 10, list = TRUE, returnTrain = T)
VD.row <- lapply(TD.row, function(x){ which( !((1:nrow(dB0t)) %in% x) )})
dB0t.TD <- lapply(TD.row, function(x){dB0t[x,]})
dG1t.TD <- lapply(TD.row, function(x){dG1t[x,]})
dB0t.VD <- lapply(VD.row, function(x){dB0t[x,]})
dG1t.VD <- lapply(VD.row, function(x){dG1t[x,]})
# TD <- list(dB0t = dB0t.td, dG1t = dG1t.td)
# VD <- list(dB0t = dB0t.vd, dG1t = dG1t.vd)
dB0t.dG1t.TD <- apply(mapply(c, dB0t.TD, dG1t.TD), 2, function(this){as.data.frame(this)})
dB0t.dG1t.VD <- apply(mapply(c, dB0t.VD, dG1t.VD), 2, function(this){as.data.frame(this)})

# 降維
Bscore.mds.TD <- lapply(dB0t.TD, function(this){
  mds1 <- vegan::metaMDS(this, distance="euclidean", k=1, trymax = 50)
  Bscore.mds <- -scores(mds1)
  return(Bscore.mds)
})
Bscore.pca.TD <- lapply(dB0t.TD, function(this){
  pca1 <- princomp(this, cor=F)
  Bscore.pca <- -pca1$scores
  return(Bscore.pca)
})
Bscore.rda.TD <- lapply(dB0t.dG1t.TD, function(this){
  this.dv <- this[, 1:4]
  this.iv <- this[, -(1:4)]
  rda1 <- vegan::rda( this.dv ~ . , data=this.iv, scale=F )
  Bscore.rda <- -summary(rda1)$sites
  return(Bscore.rda)
})
# Bscore.cca.TD <- lapply(dB0t.dG1t.TD, function(this){
#   cc1 <- candisc::cancor(cbind(HsWood, HsGrass, HsSpider, HsInsect) ~ Ga + Gb + Gc + Gd +
#       Ge + Gf + Gi + Gj + I(Ga^2) + I(Gb^2) + I(Gc^2) + I(Gd^2) +
#       I(Ge^2) + I(Gf^2) + I(Gi^2) + I(Gj^2) + (Ga + Gb + Gc + Gd +
#       Ge + Gf + Gi + Gj)^2 -Gc:Gf,
#       data=this)
#   Bscore.cca <- cc1$scores$X
#   return(Bscore.cca)
# })
Bscore.fa.TD <- lapply(dB0t.TD, function(this){
  fa1 <- factanal(this, 1, rotation="varimax", scores="regression", trace=T)
  Bscore.fa <- data.frame(FA1 = fa1$scores)
  return(Bscore.fa)
})

## multiple regression + AIC 求回歸式並以 dG1t.VD 求 prediction
# mds
lm.mds.aic.error.CV <- list()
for(i in 1:10){
  # TD 複迴歸
  # lm.cca.aic <- stats::step(lm(Bscore.mds.TD[[i]] ~ . , data=dG1t.TD[[i]]))
  # lm.cca.bic <- stats::step(lm(Bscore.mds.TD[[i]] ~ . , data=dG1t.TD[[i]]), k = log(nrow(dG1t.TD[[i]])))
  lm.mds.aic <- stats::step(lm(Bscore.mds.TD[[i]]     ~ . , data=dG1t.TD[[i]]))
  lm.mds.bic <- stats::step(lm(Bscore.mds.TD[[i]]     ~ . , data=dG1t.TD[[i]]), k = log(nrow(dG1t.TD[[i]])))
  lm.pca.aic <- stats::step(lm(Bscore.pca.TD[[i]][,1] ~ . , data=dG1t.TD[[i]]))
  lm.pca.bic <- stats::step(lm(Bscore.pca.TD[[i]][,1] ~ . , data=dG1t.TD[[i]]), k = log(nrow(dG1t.TD[[i]])))
  lm.rda.aic <- stats::step(lm(Bscore.rda.TD[[i]][,1] ~ . , data=dG1t.TD[[i]]))
  lm.rda.bic <- stats::step(lm(Bscore.rda.TD[[i]][,1] ~ . , data=dG1t.TD[[i]]), k = log(nrow(dG1t.TD[[i]])))
  lm.fa.aic  <- stats::step(lm(Bscore.fa.TD[[i]][[1]] ~ . , data=dG1t.TD[[i]]))
  lm.fa.bic  <- stats::step(lm(Bscore.fa.TD[[i]][[1]] ~ . , data=dG1t.TD[[i]]), k = log(nrow(dG1t.TD[[i]])))
  # VD 降維後分數
  score.mds.VD <- as.vector(-scores(vegan::metaMDS(dB0t.VD[[i]], distance="euclidean", k=1, trymax = 50)))
  score.pca.VD <- as.vector(-princomp(dB0t.VD[[i]], cor=F)$scores[,1])
  score.rda.VD <- as.vector(summary(vegan::rda( dB0t.VD[[i]] ~ . , data=dG1t.VD[[i]], scale=F ))$sites[,1])
  score.fa.VD  <- as.vector(factanal(dB0t.VD[[i]], 1, rotation="varimax", scores="regression", trace=T)$scores)
  # score.rda.VD <-
  #....
  lm.mds.aic.error.TD[[i]] <- score.mds.VD - as.vector(predict(lm.mds.aic, dG1t.VD[[i]]))
  #....
  # lm.cca.aic.error.TD[[i]] <- dB0t.VD[[i]] - predict(lm.cca.aic, dG1t.VD[[i]])
}
lm.mds.aic.error.CV <- unlist(lm.mds.aic.error.CV)
