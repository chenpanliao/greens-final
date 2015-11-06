y1 <- 1:11
y2 <- 1:11 + rnorm(11,0,0.1)
y3 <- 2*(1:11) + rnorm(11,0,0.1)
y4 <- runif(11)
dd <- data.frame(y1,y2,y3)

fit.pca <- princomp(dd, cor=F)
# scale(as.matrix(dd), scale=F) %*% as.matrix(loadings(fit.pca)) = fit.pca$scores 回求
# PCA分數 = 中心化原始變數 %*% loadings

fit.fa <- factanal(dd, 3, rotation="varimax", scores="regression", trace=T)
loading <- loadings(fit.fa)
score <- as.vector(fit.fa$scores)
# as.matrix(score) %*% t(as.matrix(loading)) ~= scale(dd) 回求
# FA 分數(又叫因素) %*% loadings = 標準化原始變數
# solve() solves the equation a %*% x = b for x, where b can be either a vector or a matrix.
# solve(a, b, ...)
# as.matrix(loading) %*% t(as.matrix(score)) = t(scale(dd))
# 求 FA 分數
solve(as.matrix(loading), t(scale(dd)))
solve( matrix(c(0.99,0.99,0.99,0.01,0.01,0.02,0.01,0.01,0.01), 3) , t(scale(dd)))




fit.mds <- vegan::monoMDS(dist(dd), distance="euclidean", k = 1, trymax = 21)
dd;plot(fit.mds)
score <- scores(fit.mds)
loading <- fit.mds$species
scale(scale(as.matrix(dd), center=T, scale=F) %*% as.matrix(loading)) # 回求
scale(score)
# MDS 是找到一個方向使不相似性最大
# 不相似性大概可以從 decostand(dd, "total") 看出








############# 做出所有可能假資料
areaL <- 9 # 分 areaL 級面積等級
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
par(cex=8/12)
hv <- heatmap(
  x,
  Rowv = NA,
  Colv = NA,
  col = rainbow(256),
  # col = gray(seq(1, 0, length=11)),
  scale = "column",
  # scale = "none",
  na.rm = F,
  margins = c(5,12), cexRow=1, cexCol=1,
  labRow = "",
  # labCol = c("未定義", LETTERS[1:13]),
  xlab = "", ylab =  "")
# quartz.save("./slide/所有情況模擬RdaAic.pdf", type="pdf")
quartz.save("~/Desktop/所有情況模擬RdaAic.png", dpi=300, type="png")
# dev.off()
