library(Decomp2d)

### example  : composite of two components having different frequencies
nr <- nc <- 128; x <- seq(0, 1, length=nr); y <- seq(0, 1, length=nc)

coscomp1 <- outer(cos(20 * pi * x), cos(20 * pi * y))
coscomp2 <- outer(cos(5* pi * x), cos(5 * pi * y))
cosmeanf <- coscomp1 + coscomp2
zcos0 <- range(cosmeanf)

quartz(width=6, height=2.25)
par(mfrow=c(1,3), mar=rep(0.1, 4), oma=c(0,0,1.35,0))
image(cosmeanf, xlab="", ylab="", main="", col=gray(0:100/100), axes=FALSE, zlim=zcos0)
mtext("Cosines", side=3, line=0.1, cex=0.85, font=2)
image(coscomp1, xlab="", ylab="", main="", col=gray(0:100/100), axes=FALSE, zlim=zcos0)
mtext("high-frequency component", side=3, line=0.1, cex=0.85, font=2)
image(coscomp2, xlab="", ylab="", main="", col=gray(0:100/100), axes=FALSE, zlim=zcos0)
mtext("low-frequency component", side=3, line=0.1, cex=0.85, font=2)

### decomposition by wavelet transform

# decomposition with high-frequency components above each resolution level of 2, 3, 4, 5
outcoswr2 <- decomp2d(cosmeanf, method="wavelet", wavelet.highlevel=2) 
outcoswr3 <- decomp2d(cosmeanf, method="wavelet", wavelet.highlevel=3) 
outcoswr4 <- decomp2d(cosmeanf, method="wavelet", wavelet.highlevel=4) 
outcoswr5 <- decomp2d(cosmeanf, method="wavelet", wavelet.highlevel=5) 

quartz(width=12.2/6*4, height=4.2)
par(mfcol=c(2,4), mar=rep(0.1, 4), oma=c(0,1.35,1.35,0))
image(outcoswr2$fc, xlab="", ylab="", main="",  col=gray(0:100/100), axes=FALSE, zlim=zcos0)
mtext("high-frequency component", side=2, line=0.3, cex=0.85, font=2)
mtext("level 2", side=3, line=0.1, cex=0.85, font=2)
image(outcoswr2$residue, xlab="", ylab="", main="", col=gray(0:100/100), axes=FALSE, zlim=zcos0)
mtext("low-frequency component",side=2, line=0.3, cex=0.85, font=2)
image(outcoswr3$fc, xlab="", ylab="", main="", col=gray(0:100/100), axes=FALSE, zlim=zcos0)
mtext("level 3", side=3, line=0.1, cex=0.85, font=2)
image(outcoswr3$residue, xlab="", ylab="", main="", col=gray(0:100/100), axes=FALSE, zlim=zcos0)
image(outcoswr4$fc, xlab="", ylab="", main="", col=gray(0:100/100), axes=FALSE, zlim=zcos0)
mtext("level 4", side=3, line=0.1, cex=0.85, font=2)
image(outcoswr4$residue, xlab="", ylab="", main="", col=gray(0:100/100), axes=FALSE, zlim=zcos0)
image(outcoswr5$fc, xlab="", ylab="", main="", col=gray(0:100/100), axes=FALSE, zlim=zcos0)
mtext("level 5", side=3, line=0.1, cex=0.85, font=2)
image(outcoswr5$residue, xlab="", ylab="", main="", col=gray(0:100/100), axes=FALSE, zlim=zcos0)

### decomposition by discrete cosine transform

# find most largest two coefficients of discrete cosine transform
cosmeanfdct <- DCT2D(cosmeanf, returnmat=TRUE)
coeff1 <- arrayInd(which.max(cosmeanfdct), dim(cosmeanf))
cosmeanfdct[coeff1]

cosmeanfdct[coeff1] <- 0 
coeff2 <- arrayInd(which.max(cosmeanfdct), dim(cosmeanf))  
cosmeanfdct[coeff2]

# decomposition by DCT with frequency
outcosdct <- decomp2d(cosmeanf, method="dct", dct.frequency=18)

quartz(width=6, height=2.25)
par(mfrow=c(1,3), mar=rep(0.1, 4), oma=c(0,0,1.35,0))
image(1:nrow(cosmeanf), 1:ncol(cosmeanf), abs(cosmeanfdct), xlab="", ylab="", main="", col=gray(40:100/100), axes=FALSE)
points(coeff1, col="white", pch=20, cex=1)
points(coeff2, col="white", pch=20, cex=1)
mtext("DCT coefficients", side = 3, line = 0.1, cex=0.85, font=2)
image(outcosdct$fc, xlab="", ylab="", main="", col=gray(0:100/100), axes=FALSE, zlim=zcos0)
mtext("high-frequency component", side=3, line=0.1, cex=0.85, font=2)
image(outcosdct$residue, xlab="", ylab="", main="", col=gray(0:100/100), axes=FALSE, zlim=zcos0)
mtext("low-frequency component", side=3, line= 0.1, cex=0.85, font=2)

### decomposition by two-dimensional singular spectrum analysis

# two-dimensional singular spectrum analysis with window size 64*64
outssa <- ssa(cosmeanf, kind="2d-ssa", L=c(64, 64))

quartz(width=12.2/6*4, height=4.2)
plot(outssa, type="vectors")

# weighted correlation for the reconstructed components by 10 eigenvectors
quartz()
plot(wcor(outssa, groups = 1:10), main="")

# decomposition with high-frequency component of 1 eigenvector
outssa1 <- decomp2d(cosmeanf, method="ssa", ssa.L=c(64, 64), ssa.freqcomp=1)
# decomposition with high-frequency component of 1, 2 eigenvectors
outssa2 <- decomp2d(cosmeanf, method="ssa", ssa.L=c(64, 64), ssa.freqcomp=1:2)
# decomposition with high-frequency component of 1, 2, 3 eigenvectors
outssa3 <- decomp2d(cosmeanf, method="ssa", ssa.L=c(64, 64), ssa.freqcomp=1:3)
# decomposition with high-frequency component of 1, 2, 3, 8 eigenvectors
outssa4 <- decomp2d(cosmeanf, method="ssa", ssa.L=c(64, 64), ssa.freqcomp=c(1:3, 8))

quartz(width=12.2/6*4, height=4.2)
par(mfcol=c(2,4), mar=rep(0.1, 4), oma=c(0,1.35,1.35,0))
image(outssa1$fc, xlab="", ylab="", main="", col=gray(0:100/100), axes=FALSE, zlim=zcos0)
mtext("high-frequency component",  side=2, line=0.3, cex=0.85, font=2)
mtext("1, (2, 3, 4, 5, 6, 7, 8, 9, 10)", side=3, line=0.1, cex=0.85, font=2)
image(outssa1$residue, xlab="", ylab="", main="", col=gray(0:100/100), axes=FALSE, zlim=zcos0)
mtext("low-frequency component", side=2, line=0.3, cex=0.85, font=2)
image(outssa2$fc, xlab="", ylab="", main="", col=gray(0:100/100), axes=FALSE, zlim=zcos0)
mtext("(1, 2), (3, 4, 5, 6, 7, 8, 9, 10)", side=3, line=0.1, cex=0.85, font=2)
image(outssa2$residue, xlab="", ylab="", main="", col=gray(0:100/100), axes=FALSE, zlim=zcos0)
image(outssa3$fc, xlab="", ylab="", main="", col=gray(0:100/100), axes=FALSE, zlim=zcos0)
mtext("(1, 2, 3), (4, 5, 6, 7, 8, 9, 10)", side=3, line=0.1, cex=0.85, font=2)
image(outssa3$residue, xlab="", ylab="", main="", col=gray(0:100/100), axes=FALSE, zlim=zcos0)
image(outssa4$fc, xlab="", ylab="", main="", col=gray(0:100/100), axes=FALSE, zlim=zcos0)
mtext("(1, 2, 3, 8), (4, 5, 6, 7, 9, 10)", side=3, line=0.1, cex=0.85, font=2)
image(outssa4$residue, xlab="", ylab="", main="", col=gray(0:100/100), axes=FALSE, zlim=zcos0)

### decomposition by two-dimensional principal component analysis

# contribution of principal component 
outcossvd <- svd(t(cosmeanf-mean(cosmeanf)) %*% (cosmeanf-mean(cosmeanf)))
cumsum(outcossvd$d[1:10]) / sum(outcossvd$d)

# decomposition with high-frequency component of 2 to 10 eigenvector
outcossvd1 <- decomp2d(cosmeanf, method="pca", pca.freqcomp=2:10)

quartz(width=6, height=2.25)
par(mfrow=c(1,3), mar=rep(0.1, 4), oma=c(0,2.0,1.35,0))
plot(outcossvd$d[1:10]/ sum(outcossvd$d), xlab="", ylab="", main="")
mtext("contribution", side=3, line=0.1, cex=0.85, font=2)
image(outcossvd1$fc, xlab="", ylab="", main="", col=gray(0:100/100), axes=FALSE, zlim=zcos0)
mtext("high-frequency component", side=3, line=0.1, cex=0.85, font=2)
image(outcossvd1$residue, xlab="", ylab="", main="", col=gray(0:100/100), axes=FALSE, zlim=zcos0)
mtext("low-frequency component", side=3, line=0.1, cex=0.85, font=2)

### decomposition by bidimensional empirical mode decomposition

# decomposition by conventional bidimensional empirical mode decomposition
outcosemd <- decomp2d(cosmeanf, method="emd", emd.sm=FALSE)

# decomposition by smoothing bidimensional empirical mode decomposition
outcosemd2 <- decomp2d(cosmeanf, method="emd", emd.sm=TRUE, emd.spar=0.015) 

quartz(width=12.2/6*2, height=4.2)
par(mfcol=c(2,2), mar=rep(0.1, 4), oma=c(0,1.35,1.35,0))
image(outcosemd$fc, xlab="", ylab="", main="", col=gray(0:100/100), axes=FALSE, zlim=zcos0)
mtext("high-frequency component", side=2, line=0.3, cex=0.85, font=2)
mtext("BEMD", side=3, line=0.1, cex=0.85, font=2)
image(outcosemd$residue, xlab="", ylab="", main="", col=gray(0:100/100), axes=FALSE, zlim=zcos0)
mtext("low-frequency component", side=2, line=0.3, cex=0.85, font=2)
image(outcosemd2$fc, xlab="", ylab="", main="", col=gray(0:100/100), axes=FALSE, zlim=zcos0)
mtext("BSEMD", side=3, line=0.1, cex=0.85, font=2)
image(outcosemd2$residue, xlab="", ylab="", main="", col=gray(0:100/100), axes=FALSE, zlim=zcos0)

### decomposition by ensemble patch transform

# empirical period defined by the distance of local extrema
quartz(width=9.5, height=3.5)
par(mfrow=c(1,2), mar=c(2,2,1,1))
hist(empperiod(cosmeanf)$rowperiod, xaxt = "n", breaks=seq(4, 55, by=3), freq=FALSE, 
    xlab="", main="empirical period of vertical direction")
axis(1, seq(4, 55, by=3), seq(4, 55, by=3))
hist(empperiod(cosmeanf)$colperiod, xaxt = "n", breaks=seq(4, 55, by=3), freq=FALSE, 
    xlab="", main="empirical period of horizontal direction")
axis(1, seq(4, 55, by=3), seq(4, 55, by=3))

# decomposition by ept with patch size tau = c(3,3), c(7,7), c(13,13) and c(30,30)
outcossift3 <- decomp2d(cosmeanf, method="ept", ept.tau=3)
outcossift7 <- decomp2d(cosmeanf, method="ept", ept.tau=7)
outcossift13 <- decomp2d(cosmeanf, method="ept", ept.tau=13)
outcossift30 <- decomp2d(cosmeanf, method="ept", ept.tau=30)

quartz(width=12.2/6*4, height=4.2)
par(mfcol=c(2,4), mar=rep(0.1, 4), oma=c(0,1.35,1.35,0))
image(outcossift3$fc, xlab="", ylab="", main="", col=gray(0:100/100), axes=FALSE, zlim=zcos0)
mtext("high-frequency component", side=2, line=0.3, cex=0.85, font=2)
mtext("(3, 3)", side=3, line=0.1, cex=0.85, font=2)
image(outcossift3$residue, xlab="", ylab="", main="", col=gray(0:100/100), axes=FALSE, zlim=zcos0)
mtext("low-frequency component", side=2, line=0.3, cex=0.85, font=2)
image(outcossift7$fc, xlab="", ylab="", main="", col=gray(0:100/100), axes=FALSE, zlim=zcos0)
mtext("(7, 7)", side=3, line=0.1, cex=0.85, font=2)
image(outcossift7$residue, xlab="", ylab="", main="", col=gray(0:100/100), axes=FALSE, zlim=zcos0)
image(outcossift13$fc, xlab="", ylab="", main="", col=gray(0:100/100), axes=FALSE, zlim=zcos0)
mtext("(13, 13)", side=3, line=0.1, cex=0.85, font=2)
image(outcossift13$residue, xlab="", ylab="", main="", col=gray(0:100/100), axes=FALSE, zlim=zcos0)
image(outcossift30$fc, xlab="", ylab="", main="", col=gray(0:100/100), axes=FALSE, zlim=zcos0)
mtext("(30, 30)", side=3, line=0.1, cex=0.85, font=2)
image(outcossift30$residue, xlab="", ylab="", main="", col=gray(0:100/100), axes=FALSE, zlim=zcos0)

