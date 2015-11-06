library(PMCMR)
library(caret)
library(relaimpo)
library(QuantPsyc)
library(vegan)
library(xtable)
library(magrittr)
library(tidyr)
library(ggplot2)
library(hydroGOF)

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

d <- read.csv("dat/All_Data_20151020.csv")
names(d)



# 所有樣本
# 有 NA
grep("^G[a-z]$", names(d))
dG0 <- d[, grep("^G[a-x]$", names(d))]
dG1 <- d[, grep("^G[a-w]$", names(d))]
dB0 <- d[, grep("^(Dgs|Hs)[a-zA-Z]+$", names(d))]


## 找出總面積不正常的點
quartz(width=7, height=3)
par(cex=10/12, mar=c(4,4,0,0)+0.1, family = "Noto Sans CJK TC")
hist(rowSums(dG0), nclass=40, main="", xlim=c(0,60000), ylim=c(0,500), col=3, xlab="總面積（m²）", ylab="次數")
  abline(v=c(8000,12000), lty=2, col=2)
  text(4000, 260, "< 8000\n(219/1169)")
  text(16000, 260, "> 12000\n(26/1169)")
quartz.save("./slide/invalid-area.pdf", type="pdf")
quartz.save("./slide/invalid-area.png", dpi=300, type="png")
dev.off()


## 找出未定義面積不正常的點
quartz(width=7, height=3)
par(cex=10/12, mar=c(4,4,0,0)+0.1, family = "Noto Sans CJK TC")
hist(dG0$Gx, nclass=40, main="", xlim=c(0,20000), ylim=c(0,600), col=3, xlab="未定義面積（m²）", ylab="次數")
  abline(v=c(8000), lty=2, col=2)
  # text(4000, 260, "< 8000\n(18.7%)")
  text(12000, 550, "> 8000\n(125/1169)")
quartz.save("./slide/invalid-gx-area.pdf", type="pdf")
quartz.save("./slide/invalid-gx-area.png", dpi=300, type="png")
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
  col = gray(seq(1, 0, -0.25)),
  scale = "row", na.rm = F,
  margins = c(6,6),cexRow=0.3, cexCol=1,
  labRow = Location,
  labCol = c("未定義", LETTERS[1:13]),
  xlab = "類型", ylab =  "樣點")
quartz.save("./slide/invalid-type.pdf", type="pdf")
quartz.save("./slide/invalid-type.png", dpi=300, type="png")
dev.off()


## 有問題的多樣性資料
x  <- as.matrix(dB0)
Location <- as.character(d$Location)
Location[duplicated(Location)] <- ""
quartz(width=7, height=7)
par(family = "Noto Sans CJK TC", cex=10/12)
hv <- heatmap(x,
  Rowv = NA, Colv = NA,
  # col = terrain.colors(50),
  col = gray(seq(1, 0, -0.25)),
  scale = "column", na.rm = F,
  margins = c(6,6),cexRow=0.3,cexCol=1,
  labRow = Location,
  labCol = paste(c("木本","草本","蜘蛛","昆蟲"), c(rep("Simpson",4), rep("Shannon",4)), sep="\n"),
  xlab = "生物多樣性指數", ylab =  "樣點")
quartz.save("./slide/invalid-bio.pdf", type="pdf")
quartz.save("./slide/invalid-bio.png", dpi=300, type="png")
dev.off()















# 求解用的 dataset
# 無 NA
# 多樣性資料只看 shannon
# 去除 Gl
# 去除 < 7000 | > 13000 area
# 去除 Gx > 7000
# 面積以列標準化為 sum = 1
# NO=="640" 很怪
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
# 減去沒用的交互作用
dG1t <- dG1t[!colnames(dG1t) %in% c("GaGc","GaGf","GbGc","GbGf","GcGd","GcGe","GcGf","GcGi","GcGj","GdGf","GeGf","GfGi","GfGj","Gc2","Gf2")]
dG1t.nointer <- dG1t[, 1:8]
dG1t.order2  <- dG1t[, 24:29]

# 生態資料
dB0t <- decostand(dt[, c(22,23,24,25)], "range", 2)
  # dB0t <- decostand(dt[, c(22,23,24,25)], 2, "standardize")
names(dB0t); dim(dB0t) # 783 spots remained

# 所有採用資料圖
quartz(width=7, height=3.5)
par(cex=8/12, mar=c(5,5,0,0)+0.1, family = "Noto Sans CJK TC", mgp=c(4,1,0))
boxplot(cbind(dB0t, dG1t), xlab="變數", ylab="標準化後數值", las=2)
  # axis(1, bp, names(dG1t))
quartz.save("./slide/var-summary.pdf", type="pdf")
quartz.save("./slide/var-summary.png", dpi=300, type="png")
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
  fit <- vegan::metaMDS(dt, distance="euclidean", k = 1, trymax = 20)
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
quartz.save("./slide/原始資料降維分數.png", dpi=300, type="png")
quartz.save("./slide/原始資料降維分數.pdf", type="pdf")
dev.off()







## multiple regression: diversity ~ area^1 + area^2 + interaction(among area^1)
## AIC
# dG1t.nointer <- dG1t[, 1:8]
# dG1t.order2  <- dG1t[, 24:29]
# IVsFormula <- paste(c(paste(c(
#   names(as.data.frame(model.matrix(~(.)^2, dG1t.nointer[1,])))[-1],
#   names(as.data.frame(
#     # model.matrix(~., dG1t.nointer^2[1,])
#   ))[-1]
#   ), sep="", collapse="+") ,
#   c("-GaGc","-GaGf","-GbGc","-GbGf","-GcGd","-GcGe","-GcGf","-GcGi","-GcGj", "-GdGf","-GeGf","-GfGi","-GfGj","-Gc2","-Gf2")
#   ), sep="", collapse="")
# lm(Bscore.o$PCA ~ as.formula(IVsFormula) , data=dG1t)
fit.o <- list()
fit.o$PCA.AIC <- stats::step(lm(Bscore.o$PCA ~ . , data=dG1t), scope=list(upper = ~., lower = ~1))
fit.o$PCA.BIC <- stats::step(lm(Bscore.o$PCA ~ . , data=dG1t), scope=list(upper = ~., lower = ~1), k = log(nrow(dG1t)))
fit.o$FA.AIC  <- stats::step(lm(Bscore.o$FA  ~ . , data=dG1t), scope=list(upper = ~., lower = ~1))
fit.o$FA.BIC  <- stats::step(lm(Bscore.o$FA  ~ . , data=dG1t), scope=list(upper = ~., lower = ~1), k = log(nrow(dG1t)))
fit.o$MDS.AIC <- stats::step(lm(Bscore.o$MDS ~ . , data=dG1t), scope=list(upper = ~., lower = ~1))
fit.o$MDS.BIC <- stats::step(lm(Bscore.o$MDS ~ . , data=dG1t), scope=list(upper = ~., lower = ~1), k = log(nrow(dG1t)))
fit.o$RDA.AIC <- stats::step(lm(Bscore.o$RDA ~ . , data=dG1t), scope=list(upper = ~., lower = ~1))
fit.o$RDA.BIC <- stats::step(lm(Bscore.o$RDA ~ . , data=dG1t), scope=list(upper = ~., lower = ~1), k = log(nrow(dG1t)))
#   tmp <-
#   ordiR2step(
#     vegan::rda(dB0t ~ 1 , data=dG1t, scale=F),
#     formula(vegan::rda(dB0t ~ . , data=dG1t, scale=F)),
#     Pin=0.10, Pout=0.15, permutations = how(nperm = 499)
#   )
#   paste(attr(formula(tmp), "term.labels"), sep="", collapse="+")
# fit.o$RDA.R2a <- lm(tmp$CCA$wa[,1] ~ Ge+Gd+Gj+GaGd+Ga+Gd2+GeGi+GdGi+GeGj+Ga2+Gf+Ge2+GaGi+Gc+Gb+GbGd+GaGe+Gi, data = dG1t)

# 結合係數成data.frame: coefM.o
{coefM.o <- multimerge(lapply(fit.o, function(x){data.frame(coef(x))}))} %>% print %>% xtable(., digits=4)
  # coefM.o <- merge(coefM.o, as.data.frame(coef(lm(dB0t$HsGrass ~ . , data=dG1t))), by="row.names", sort=T, all=T)
  # rownames(coefM.o) <- coefM.o$Row.names
  # coefM.o$`Row.names` <- NULL
  # coefM.o[, ncol(coefM.o)] <- NULL
coefM.o %>% print %>% xtable
rbind(
  `# IV` = sapply(fit.o, FUN = function(fit) {length(attr(fit$terms, "term.labels"))} ),
  NRMSE = sapply(fit.o, FUN = function(fit) {sqrt(sum(resid(fit)^2) / length(fit$residuals)) / sd(fit$model[[1]])} ),
  R2 = sapply(fit.o, FUN = function(fit) {summary(fit)$r.squared} ),
  `Adjusted R2` = sapply(fit.o, FUN = function(fit) {summary(fit)$adj.r.squared} ),
  `AIC/1000` = sapply(fit.o, AIC) / 1000,
  `BIC/1000` = sapply(fit.o, BIC) / 1000
) %>% print %>% xtable(., digits=3)

quartz(width=7, height=2)
par(mar=c(4,4,0,0)+0.1, family = "Noto Sans CJK TC", mgp=c(3,1,0), cex=8/12)
boxplot(sapply(fit.o, resid), xlab="模型", ylab="殘差", las=1)
quartz.save("./slide/resid最佳迴歸法.png", dpi=300, type="png")
quartz.save("./slide/resid最佳迴歸法.pdf", type="pdf")
dev.off()
# quartz(width=7, height=2)
# tmp <- fit.o %>% sapply(., resid) %>% as.data.frame %>% gather(., Model, Residual, 1:8)
# ggplot(tmp, aes(Model, Residual)) + geom_jitter(colour=gray(0.5),size=1,height=10) + geom_violin(alpha=0.8) + theme_bw()


# 結合標準化係數成data.frame: coefM.o.std 使用 QuantPsyc::lm.beta()
# coefM.o.std <- multimerge(lapply(fit.o, function(x){data.frame(lm.beta(x))}))
#   coefM.o.std <- merge(coefM.o.std, as.data.frame(coef(lm(dB0t$HsGrass ~ . , data=dG1t))), by="row.names", sort=T, all=T)
#   rownames(coefM.o.std) <- coefM.o.std$Row.names
#   coefM.o.std$`Row.names` <- NULL
#   coefM.o.std[, ncol(coefM.o.std)] <- NULL
# round(coefM.o.std, digit=3)
# xtable(coefM.o.std)



## beta : relative importance
ri <- list()
for(i in 1:length(fit.o)) {
  ri[[i]] <- calc.relimp(fit.o[[i]], type=c("lmg"), rela=F)
}
quartz(width=10, height=5)
par(mar=c(5,5,1,0)+0.1, family = "Noto Sans CJK TC", mgp=c(4,1,0), cex=4/12, mfrow=c(2,4))
for(i in 1:length(ri)) {
  barplot( ri[[i]]@lmg, xlab="", ylab="R² 貢獻量", beside=T, ylim=c(0,0.2), las=2)
  title(names(fit.o)[i])
}
quartz.save( "beta相對貢獻量.pdf", type="pdf")
quartz.save( "beta相對貢獻量.png", type="png", dpi=300)
dev.off()
# boot <- boot.relimp(fit.o[[1]], b = 1000, type = c("lmg"), rank = TRUE, diff = TRUE, rela = TRUE)
# booteval.relimp(boot) # print result
# plot(booteval.relimp(boot,sort=F)) # plot result




## beta plot before standardizing
tmp <- data.frame(
  aiccol = 1:4,
  biccol = 5:8,
  filename = paste0("beta-", sub(".AIC", "", colnames(coefM.o)[1:4]))
)
for(i in 1:nrow(tmp)) {
  quartz(width=7, height=3.5)
  par(mar=c(5,5,0,0)+0.1, family = "Noto Sans CJK TC", mgp=c(4,1,0), cex=8/12)
  plot(NULL, xlim=c(1, nrow(coefM.o)), ylim=c(-15,15), xaxt="n", yaxt="n", xlab="類型", ylab="原始迴歸係數" )
    axis(1, 1:nrow(coefM.o), rownames(coefM.o), las=2)
    axis(2, seq(-16, 16, 2), seq(-16, 16, 2), las=2)
    abline(h=0)
    arrows(
      1:nrow(coefM.o)-0.2,
      rep(0, nrow(coefM.o)),
      1:nrow(coefM.o)-0.2,
      ifelse(
        coefM.o[, tmp$aiccol[i]] >= 0 | is.na(coefM.o[, tmp$aiccol[i]]),
        coefM.o[,tmp$aiccol[i]],
        coefM.o[,tmp$aiccol[i]]
      ),
      code=2, length = 0.05, col=1 #, lwd=abs(coefM.o[,1])/max(abs(coefM.o[,1]))*5+1
    )
    arrows(
      1:nrow(coefM.o)+0.2,
      rep(0, nrow(coefM.o)),
      1:nrow(coefM.o)+0.2,
      ifelse(
        coefM.o[,tmp$biccol[i]] >= 0 | is.na(coefM.o[,tmp$biccol[i]]),
        coefM.o[,tmp$biccol[i]],
        coefM.o[,tmp$biccol[i]]
      ) ,
      code=2, length = 0.05, col=2, lty=1#, lwd=abs(coefM.o[,2])/max(abs(coefM.o[,2]))*5+1
    )
    legend(1, 8, c("AIC", "BIC"), lty=c(1,1), col=c(1,2))
  quartz.save( paste0(tmp$filename[i], ".png") , dpi=300, type="png")
  quartz.save( paste0(tmp$filename[i], ".pdf"), type="pdf")
  dev.off()
}


## beta plot after standardizing
tmp <- data.frame(
  aiccol = 1:4,
  biccol = 5:8,
  filename = paste0("beta-std", sub(".AIC", "", colnames(coefM.o.std)[1:4]))
)
for(i in 1:nrow(tmp)) {
  quartz(width=7, height=3.5)
  par(mar=c(5,5,0,0)+0.1, family = "Noto Sans CJK TC", mgp=c(4,1,0), cex=8/12)
  plot(NULL, xlim=c(1, nrow(coefM.o.std)), ylim=c(-0.6,1), xaxt="n", yaxt="n", xlab="類型", ylab="標準化迴歸係數" )
    axis(1, 1:nrow(coefM.o.std), rownames(coefM.o.std), las=2)
    axis(2, seq(-0.6, 1, 0.2), round(seq(-0.6, 1, 0.2),2), las=2)
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
  quartz.save( paste0(tmp$filename[i], ".png") , dpi=300, type="png")
  quartz.save( paste0(tmp$filename[i], ".pdf"), type="pdf")
  dev.off()
}




########### find best fitting BAF
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
quartz.save("./slide/all.png", dpi=300, type="png")
quartz.save("./slide/all.pdf", type="pdf")
dev.off()















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
  fit.TD[[i]]$PCA.AIC <- stats::step(lm(Bscore.TD[[i]]$PCA ~ . , data=dG1t.TD[[i]]), scope=list(upper = ~., lower = ~1))
  fit.TD[[i]]$PCA.BIC <- stats::step(lm(Bscore.TD[[i]]$PCA ~ . , data=dG1t.TD[[i]]), scope=list(upper = ~., lower = ~1), k = log(nrow(dG1t.TD[[i]])))
  fit.TD[[i]]$FA.AIC  <- stats::step(lm(Bscore.TD[[i]]$FA  ~ . , data=dG1t.TD[[i]]), scope=list(upper = ~., lower = ~1))
  fit.TD[[i]]$FA.BIC  <- stats::step(lm(Bscore.TD[[i]]$FA  ~ . , data=dG1t.TD[[i]]), scope=list(upper = ~., lower = ~1), k = log(nrow(dG1t.TD[[i]])))
  fit.TD[[i]]$MDS.AIC <- stats::step(lm(Bscore.TD[[i]]$MDS ~ . , data=dG1t.TD[[i]]), scope=list(upper = ~., lower = ~1))
  fit.TD[[i]]$MDS.BIC <- stats::step(lm(Bscore.TD[[i]]$MDS ~ . , data=dG1t.TD[[i]]), scope=list(upper = ~., lower = ~1), k = log(nrow(dG1t.TD[[i]])))
  fit.TD[[i]]$RDA.AIC <- stats::step(lm(Bscore.TD[[i]]$RDA ~ . , data=dG1t.TD[[i]]), scope=list(upper = ~., lower = ~1))
  fit.TD[[i]]$RDA.BIC <- stats::step(lm(Bscore.TD[[i]]$RDA ~ . , data=dG1t.TD[[i]]), scope=list(upper = ~., lower = ~1), k = log(nrow(dG1t.TD[[i]])))

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
rbind(NRMSE = sqrt(colSums((Bscore.TD.pred.求最佳迴歸法 - Bscore.VD.求最佳迴歸法)^2) / nrow(Bscore.TD.pred.求最佳迴歸法)) / apply(Bscore.VD.求最佳迴歸法, 2, sd))  %>% print %>% xtable(., digits=4)
# seed = 52004800
# PCA.AIC.p PCA.BIC.p  FA.AIC.p  FA.BIC.p MDS.AIC.p MDS.BIC.p RDA.AIC.p RDA.BIC.p
# NRMSE 0.7537227 0.7585283 0.7866131 0.7888567 0.7584514 0.7609845 0.7378838  0.738952
# % latex table generated in R 3.2.2 by xtable 1.7-4 package
# % Fri Nov  6 01:18:50 2015
# \begin{table}[ht]
# \centering
# \begin{tabular}{rrrrrrrrr}
# \hline
# & PCA.AIC.p & PCA.BIC.p & FA.AIC.p & FA.BIC.p & MDS.AIC.p & MDS.BIC.p & RDA.AIC.p & RDA.BIC.p \\
# \hline
# NRMSE & 0.7537 & 0.7585 & 0.7866 & 0.7889 & 0.7585 & 0.7610 & 0.7379 & 0.7390 \\
# \hline
# \end{tabular}
# \end{table}
quartz(width=7, height=2.5)
par(mar=c(4,4,0,0)+0.1, family = "Noto Sans CJK TC", mgp=c(3,1,0), cex=8/12)
boxplot(  Bscore.TD.pred.求最佳迴歸法 - Bscore.VD.求最佳迴歸法 , xlab="模型", ylab="訓練與驗證差值", las=1)
quartz.save("./slide/CV求最佳迴歸法.png", dpi=300, type="png")
quartz.save("./slide/CV求最佳迴歸法.pdf", type="pdf")
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
  Bscore.TD[[i]]$PCA.AIC <- Bscore.TD[[i]]$PCA.BIC <- Bscore.o[TD.row[[i]],]$PCA
  Bscore.TD[[i]]$FA.AIC  <- Bscore.TD[[i]]$FA.BIC  <- Bscore.o[TD.row[[i]],]$FA
  Bscore.TD[[i]]$MDS.AIC <- Bscore.TD[[i]]$MDS.BIC <- Bscore.o[TD.row[[i]],]$MDS
  Bscore.TD[[i]]$RDA.AIC <- Bscore.TD[[i]]$RDA.BIC <- Bscore.o[TD.row[[i]],]$RDA
  rownames(Bscore.TD[[i]]) <- dt$NO[TD.row[[i]]]
  Bscore.TD[[i]]$NO <- NULL

  # 生態驗證資料 dB0t.VD[[i]] = Bscore.VD[[i]]          (x)
  # 已降維生態驗證資料 dB0t.VD[[i]] 命為 Bscore.VD[[i]] (o)
  Bscore.VD[[i]] <- data.frame(NO = dt$NO[VD.row[[i]]])
  Bscore.VD[[i]]$PCA.AIC <- Bscore.VD[[i]]$PCA.BIC <- Bscore.o[VD.row[[i]],]$PCA
  Bscore.VD[[i]]$FA.AIC  <- Bscore.VD[[i]]$FA.BIC  <- Bscore.o[VD.row[[i]],]$FA
  Bscore.VD[[i]]$MDS.AIC <- Bscore.VD[[i]]$MDS.BIC <- Bscore.o[VD.row[[i]],]$MDS
  Bscore.VD[[i]]$RDA.AIC <- Bscore.VD[[i]]$RDA.BIC <- Bscore.o[VD.row[[i]],]$RDA
  rownames(Bscore.VD[[i]]) <- dt$NO[VD.row[[i]]]
  Bscore.VD[[i]]$NO <- NULL

  # 拿 Bscore.TD[[i]]$xxA.XIC 於 formu.o$xxA.xIC 做複迴歸
  fit.TD[[i]] <- list()
  fit.TD[[i]]$PCA.AIC <- lm(paste( "Bscore.TD[[i]]$PCA.AIC ~" , formu.o$PCA.AIC) , data=dG1t.TD[[i]])
  fit.TD[[i]]$PCA.BIC <- lm(paste( "Bscore.TD[[i]]$PCA.BIC ~" , formu.o$PCA.BIC) , data=dG1t.TD[[i]])
  fit.TD[[i]]$FA.AIC  <- lm(paste( "Bscore.TD[[i]]$FA.AIC ~"  , formu.o$FA.AIC)  , data=dG1t.TD[[i]])
  fit.TD[[i]]$FA.BIC  <- lm(paste( "Bscore.TD[[i]]$FA.BIC ~"  , formu.o$FA.BIC)  , data=dG1t.TD[[i]])
  fit.TD[[i]]$MDS.AIC <- lm(paste( "Bscore.TD[[i]]$MDS.AIC ~" , formu.o$MDS.AIC) , data=dG1t.TD[[i]])
  fit.TD[[i]]$MDS.BIC <- lm(paste( "Bscore.TD[[i]]$MDS.BIC ~" , formu.o$MDS.BIC) , data=dG1t.TD[[i]])
  fit.TD[[i]]$RDA.AIC <- lm(paste( "Bscore.TD[[i]]$RDA.AIC ~" , formu.o$RDA.AIC) , data=dG1t.TD[[i]])
  fit.TD[[i]]$RDA.BIC <- lm(paste( "Bscore.TD[[i]]$RDA.BIC ~" , formu.o$RDA.BIC) , data=dG1t.TD[[i]])

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
Bscore.TD.pred.最佳解釋變數組 <- do.call("rbind", Bscore.TD.pred)
  Bscore.TD.pred.最佳解釋變數組 <- Bscore.TD.pred.final[order(as.numeric(rownames(Bscore.TD.pred.final))), ]
Bscore.VD.最佳解釋變數組 <- do.call("rbind", Bscore.VD)
  Bscore.VD.最佳解釋變數組 <- Bscore.VD.final[order(as.numeric(rownames(Bscore.VD.final))), ]
dim(Bscore.TD.pred.最佳解釋變數組) ; dim(Bscore.VD.最佳解釋變數組)
names(Bscore.TD.pred.最佳解釋變數組) ; names(Bscore.VD.最佳解釋變數組)
# 看看誰最好
rbind(NRMSE = sqrt(colSums((Bscore.TD.pred.最佳解釋變數組 - Bscore.VD.最佳解釋變數組)^2) / nrow(Bscore.TD.pred.最佳解釋變數組)) / apply(Bscore.VD.最佳解釋變數組, 2, sd))  %>% print %>% xtable(., digits=4)
# PCA.AIC.p PCA.BIC.p  FA.AIC.p  FA.BIC.p     NMDS1     NMDS1 RDA.AIC.p RDA.BIC.p
# NRMSE 0.7433187 0.7447574 0.8087097 0.8287309 0.7492935 0.7420321 0.7247877 0.7258387
# % latex table generated in R 3.2.2 by xtable 1.7-4 package
# % Fri Nov  6 01:46:19 2015
# \begin{table}[ht]
# \centering
# \begin{tabular}{rrrrrrrrr}
# \hline
# & PCA.AIC.p & PCA.BIC.p & FA.AIC.p & FA.BIC.p & NMDS1 & NMDS1 & RDA.AIC.p & RDA.BIC.p \\
# \hline
# NRMSE & 0.7433 & 0.7448 & 0.8087 & 0.8287 & 0.7493 & 0.7420 & 0.7248 & 0.7258 \\
# \hline
# \end{tabular}
# \end{table}
quartz(width=7, height=2.5)
par(mar=c(4,4,0,0)+0.1, family = "Noto Sans CJK TC", mgp=c(3,1,0), cex=8/12)
boxplot(  Bscore.TD.pred.最佳解釋變數組 - Bscore.VD.最佳解釋變數組 , xlab="模型", ylab="訓練與驗證差值", lax=1)
quartz.save("./slide/CV-最佳解釋變數組.png", dpi=300, type="png")
quartz.save("./slide/CV-最佳解釋變數組.pdf", type="pdf")
dev.off()











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
  fit.TD[[i]]$PCA.AIC <- stats::step(lm(Bscore.TD[[i]]$PCA ~ . , data=dG1t.TD[[i]]), scope=list(upper = ~., lower = ~1))
  fit.TD[[i]]$PCA.BIC <- stats::step(lm(Bscore.TD[[i]]$PCA ~ . , data=dG1t.TD[[i]]), scope=list(upper = ~., lower = ~1), k = log(nrow(dG1t.TD[[i]])))
  fit.TD[[i]]$FA.AIC  <- stats::step(lm(Bscore.TD[[i]]$FA  ~ . , data=dG1t.TD[[i]]), scope=list(upper = ~., lower = ~1))
  fit.TD[[i]]$FA.BIC  <- stats::step(lm(Bscore.TD[[i]]$FA  ~ . , data=dG1t.TD[[i]]), scope=list(upper = ~., lower = ~1), k = log(nrow(dG1t.TD[[i]])))
  fit.TD[[i]]$MDS.AIC <- stats::step(lm(Bscore.TD[[i]]$MDS ~ . , data=dG1t.TD[[i]]), scope=list(upper = ~., lower = ~1))
  fit.TD[[i]]$MDS.BIC <- stats::step(lm(Bscore.TD[[i]]$MDS ~ . , data=dG1t.TD[[i]]), scope=list(upper = ~., lower = ~1), k = log(nrow(dG1t.TD[[i]])))
  fit.TD[[i]]$RDA.AIC <- stats::step(lm(Bscore.TD[[i]]$RDA ~ . , data=dG1t.TD[[i]]), scope=list(upper = ~., lower = ~1))
  fit.TD[[i]]$RDA.BIC <- stats::step(lm(Bscore.TD[[i]]$RDA ~ . , data=dG1t.TD[[i]]), scope=list(upper = ~., lower = ~1), k = log(nrow(dG1t.TD[[i]])))

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
rbind(NRMSE = sqrt(colSums((Bscore.TD.pred.final - Bscore.VD.final)^2) / nrow(Bscore.TD.pred.final)) / apply(Bscore.VD.final, 2, sd))  %>% print %>% xtable(., digits=4)
# seed = 52004800
# PCA.AIC.p PCA.BIC.p  FA.AIC.p  FA.BIC.p    NMDS1     NMDS1 RDA.AIC.p RDA.BIC.p
# NRMSE 0.7599851 0.7670764 0.8742578 0.8773687 0.854141 0.8574579 0.7588922  0.760086
# % latex table generated in R 3.2.2 by xtable 1.7-4 package
# % Fri Nov  6 01:57:06 2015
# \begin{table}[ht]
# \centering
# \begin{tabular}{rrrrrrrrr}
# \hline
# & PCA.AIC.p & PCA.BIC.p & FA.AIC.p & FA.BIC.p & NMDS1 & NMDS1 & RDA.AIC.p & RDA.BIC.p \\
# \hline
# NRMSE & 0.7600 & 0.7671 & 0.8743 & 0.8774 & 0.8541 & 0.8575 & 0.7589 & 0.7601 \\
# \hline
# \end{tabular}
# \end{table}
quartz(width=7, height=2.5)
par(mar=c(4,4,0,0)+0.1, family = "Noto Sans CJK TC", mgp=c(3,1,0), cex=8/12)
boxplot(  Bscore.TD.pred.final - Bscore.VD.final , xlab="模型", ylab="訓練與驗證差值", lax=1)
quartz.save("./slide/CV-final.png", dpi=300, type="png")
quartz.save("./slide/CV-final.pdf", type="pdf")
dev.off()
# max(Bscore.TD.pred.final[,1])
# dt [max(Bscore.TD.pred.final[,1]) == Bscore.TD.pred.final[,1] , ]
# kruskal.test(value ~ Var2,  data=melt(tmp))
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









########## 印出迴歸式
coefM.o[c("RDA.AIC", "RDA.BIC")] %>% xtable(., digits=4) %>% print(., NA.string="---")


##################### 亂數模型

length(pred.pca.aic <- predict(fit.o$PCA.AIC, dG1t))
length(pred.rda.bic <- predict(fit.o$RDA.BIC, dG1t))
dt.pred <- data.frame( dt, pred.pca.aic, pred.rda.bic  )
head(dt.pred[order(-dt.pred$pred.pca.aic), ])

#### 創造假資料
rown <- 500000
dG1t.simu0 <- dG1t[rep(1,rown),1:8]
rownames(dG1t.simu0) <- NULL
set.seed(52004800)
dG1t.simu0[,] <- runif(ncol(dG1t.simu0) * rown)
dG1t.simu0 <- decostand(dG1t.simu0, "total")
head(dG1t.simu0)

# dG1t.simu0 完成：rowSums() 為 1
dG1t.simu1 <- as.data.frame(apply(dG1t.simu0, 2, function(x){x-mean(x)}))
tmp <- as.data.frame(dG1t.simu1^2); names(tmp) <- paste0(names(tmp), "2")
dG1t.simu2 <- cbind(as.data.frame(model.matrix(~(.)^2, dG1t.simu1)), tmp)  [,-1]
names(dG1t.simu2) <- gsub(":", "", names(dG1t.simu2))
names(dG1t.simu2);colMeans(dG1t.simu2);dim(dG1t.simu2)
  # 減去沒用的交互作用 -GbGc -GbGf -GcGd -GcGf -GcGi -GcGj -GfGi -GfGj
dG1t.simu3 <- dG1t.simu2[!colnames(dG1t.simu2) %in% c("GaGc","GaGf","GbGc","GbGf","GcGd","GcGe","GcGf","GcGi","GcGj","GdGf","GeGf","GfGi","GfGj","Gc2","Gf2")]
dim(dG1t.simu3);head(dG1t.simu3)
# dG1t.simu3 完成：已建好二次式和交互作用

#### 代入迴歸式
pred.rda.aic <- predict(fit.o$RDA.AIC, dG1t.simu3)
pred.rda.bic <- predict(fit.o$RDA.BIC, dG1t.simu3)

#### Top 20 sites pca.aic
t.first <- cbind(
  Bio = pred.rda.aic[order(-pred.rda.aic)[1:20]],
  dG1t.simu0[order(-pred.rda.aic)[1:20], ]*100
)
#### Top -20 sites
t.last <- cbind(
  Bio = pred.rda.aic[order(pred.rda.aic)[1:20]],
  dG1t.simu0[order(pred.rda.aic)[1:20], ]*100
)
xtable(round(t.first, 2))
xtable(round(t.last,  2))

#### Top 20 sites rda.bic
t.first <- cbind(
  Bio = pred.rda.bic[order(-pred.rda.bic)[1:20]],
  dG1t.simu0[order(-pred.rda.bic)[1:20], ]*100
)
#### Top -20 sites
t.last <- cbind(
  Bio = pred.rda.bic[order(pred.rda.bic)[1:20]],
  dG1t.simu0[order(pred.rda.bic)[1:20], ]*100
)
xtable(round(t.first, 2))
xtable(round(t.last,  2))










############# 做出所有可能假資料
areaL <- 6 # 分 areaL 級面積等級
typeN <- 6  # 有 typeN 種面積 a b d e i j
sizeM <- areaL ^ (typeN - 1)
tmp0 <- vector("list", typeN)
names(tmp0) <- paste0("G", c("a","b","d","e","i","j"))
for(i in 1:typeN){
  tmp0[[i]] <- as.numeric(gl(areaL, areaL^(i-1), sizeM)) - 1
}
dGsimAll1 <- do.call("cbind", tmp0) %>% data.frame(., check.names = F, stringsAsFactors = F)
dGsimAll1$Gj <- sum(1:typeN) - rowSums(dGsimAll1)
dGsimAll1 %<>% subset(., Gj >= 0) %>% unique %>% decostand(., "total", 1) %>% apply(., 2, function(x){x-mean(x)}) %>% as.data.frame
head(dGsimAll1) ; rowSums(dGsimAll1)
# dGsimAll1 完成：已建好一次式
dGsimAll2 <- as.data.frame(dGsimAll1^2); names(dGsimAll2) <- paste0(names(dGsimAll2), "2")
# dGsimAll2 完成：已建好二次式和交互作用
dGsimAll <- cbind(as.data.frame(model.matrix(~(.)^2, dGsimAll1)), dGsimAll2)  [,-1]
names(dGsimAll) <- gsub(":", "", names(dGsimAll))
dGsimAll <- dGsimAll[names(fit.o$RDA.AIC$model)[-1]] # 只挑出 fit.o$RDA.AIC 用得到的IV變數
dGsimAll <- dGsimAll[, order(names(dGsimAll))]
head(dGsimAll);dim(dGsimAll)
# dGsimAll 完成：結合一次式、二次式和交互作用
#### 代入迴歸式
pred.dGsimAll.rda.aic <- predict(fit.o$RDA.AIC, dGsimAll)
pred.dGsimAll.rda.bic <- predict(fit.o$RDA.BIC, dGsimAll)

## heatmap RDA.AIC
# x <- as.matrix(data.frame(Prediction = pred.dGsimAll.rda.aic/10, dGsimAll)) [ order(pred.dGsimAll.rda.aic), ]
x <- data.frame(Prediction = pred.dGsimAll.rda.aic/10, dGsimAll) %>% as.matrix %>% .[order(-pred.dGsimAll.rda.aic),]
colnames(x) %<>% gsub("G([a-z]{1})", "\\1", .)
colnames(x)[-1] %<>% toupper
quartz(width=7, height=7)
par(family = "Noto Sans CJK TC", cex=10/12)
hv <- heatmap(
  x,
  Rowv = NA, Colv = NA,
  # col = terrain.colors(50),
  col = gray(seq(1, 0, length=3)),
  scale = "column",
  # scale = "none",
  na.rm = F,
  margins = c(5,2), cexRow=1, cexCol=1,
  labRow = "",
  # labCol = c("未定義", LETTERS[1:13]),
  xlab = "", ylab =  "預測多樣性（基於 RDA.AIC）及模擬樣點面積（越黑表示同欄內相對值越大）")
quartz.save("./slide/所有情況模擬RdaAic.pdf", type="pdf")
quartz.save("./slide/所有情況模擬RdaAic.png", dpi=300, type="png")
dev.off()

## heatmap RDA.BIC
# x <- as.matrix(data.frame(Prediction = pred.dGsimAll.rda.bic/10, dGsimAll)) [ order(pred.dGsimAll.rda.bic), ]
x <- data.frame(Prediction = pred.dGsimAll.rda.bic/10, dGsimAll) %>% decostand(., "standardize", 2) %>% as.matrix %>% .[order(-pred.dGsimAll.rda.bic),]
colnames(x) %<>% gsub("G([a-z]{1})", "\\1", .)
colnames(x)[-1] %<>% toupper
quartz(width=7, height=7)
par(family = "Noto Sans CJK TC", cex=10/12)
hv <- heatmap(
  x,
  Rowv = NA, Colv = NA,
  # col = terrain.colors(50),
  col = gray(seq(1, 0, length=3)),
  scale = "column",
  # scale = "none",
  na.rm = F,
  margins = c(5,2), cexRow=1, cexCol=1,
  labRow = "",
  # labCol = c("未定義", LETTERS[1:13]),
  xlab = "", ylab =  "預測多樣性（基於 RDA.BIC）及模擬樣點面積（越黑表示同欄內相對值越大）")
quartz.save("./slide/所有情況模擬RdaBic.pdf", type="pdf")
quartz.save("./slide/所有情況模擬RdaBic.png", dpi=300, type="png")
dev.off()
