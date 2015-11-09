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
