library(Decomp2d)

### D11 texture
data(d11); zd11 <- range(d11)
par(mar=c(0,0,0,0))
image(d11, xlab="", ylab="", main="", col=gray(0:100/100), axes=FALSE)

### decomposition by wavelet transform
# decomposition with high-frequency components above resolution level of 3
outd11wr3 <- decomp2d(d11, method="wavelet", wavelet.highlevel=3)

### decomposition by discrete cosine transform
# decomposition by DCT with frequency
outd11dct <- decomp2d(d11, method="dct", dct.frequency=15)

### decomposition by two-dimensional singular spectrum analysis
# two-dimensional singular spectrum analysis with window size 64*64
outssa <- ssa(d11, kind="2d-ssa", L=c(64, 64))

quartz(width=12.2/6*4, height=4.2)
plot(outssa, type="vectors", numvectors=20)

# weighted correlation for the reconstructed components by 20 eigenvectors
quartz()
plot(wcor(outssa, groups = 1:20), main="")

# decomposition with low-frequency component of 1, 2, 3 eigenvectors
outssa2 <- decomp2d(d11, method="ssa", ssa.L=c(64, 64), ssa.freqcomp=1:3)

### decomposition by two-dimensional principal component analysis
# contribution of principal component 
outd11svd <- svd(t(d11-mean(d11)) %*% (d11-mean(d11)))
outd11svd$d[1:20] / sum(outd11svd$d)

# decomposition with frequency component
outd11svd1 <- decomp2d(d11, method="pca", pca.freqcomp=1)

### decomposition by bidimensional empirical mode decomposition
# decomposition by smoothing bidimensional empirical mode decomposition
outd11emd2 <- decomp2d(d11, method="emd", emd.sm=TRUE, emd.spar=0.01)

### decomposition by ensemble patch transform
# decomposition by ept with patch size tau = c(14,14)
outd11sift14 <- decomp2d(d11, method="ept", ept.tau=14)

### Decomposition results
quartz(width=4.2, height=12.2/6*7)
par(mfrow=c(7,2), mar=c(0.1,0.1,0.1,0.1), oma=c(0,1.35,1.35,0))

image(d11, xlab="", ylab="", main="", col=gray(0:100/100), axes=FALSE, zlim=zd11)
mtext("high-frequency component", side = 3, line = 0.1, cex=0.85, font=2)
mtext("D11", side=2, line=0.3, cex=0.85, font=2)
plot(1, type="n", axes=FALSE)
mtext("low-frequency component", side = 3, line = 0.1, cex=0.85, font=2)

image(outd11wr3$fc, xlab="", ylab="", main="", col=gray(0:100/100), axes=FALSE, zlim=zd11)
mtext("Wavelet", side=2, line=0.3, cex=0.85, font=2)
image(outd11wr3$residue, xlab="", ylab="", main="", col=gray(0:100/100), axes=FALSE, zlim=zd11)

image(outd11dct$fc, xlab="", ylab="", main="", col=gray(0:100/100), axes=FALSE, zlim=zd11)
mtext("DCT", side=2, line=0.3, cex=0.85, font=2)
image(outd11dct$residue, xlab="", ylab="", main="", col=gray(0:100/100), axes=FALSE, zlim=zd11)

image(outssa2$residue, xlab="", ylab="", main="", col=gray(0:100/100), axes=FALSE, zlim=zd11)
mtext("2DSSA", side=2, line=0.3, cex=0.85, font=2)
image(outssa2$fc, xlab="", ylab="", main="", col=gray(0:100/100), axes=FALSE, zlim=zd11)

image(outd11svd1$residue, xlab="", ylab="", main="", col=gray(0:100/100), axes=FALSE, zlim=zd11)
mtext("2DPCA", side=2, line=0.3, cex=0.85, font=2)
image(outd11svd1$fc, xlab="", ylab="", main="", col=gray(0:100/100), axes=FALSE, zlim=zd11)

image(outd11emd2$fc, xlab="", ylab="", main="", col=gray(0:100/100), axes=FALSE, zlim=zd11)
mtext("BSEMD", side=2, line=0.3, cex=0.85, font=2)
image(outd11emd2$residue, xlab="", ylab="", main="",  col=gray(0:100/100), axes=FALSE, zlim=zd11)

image(outd11sift14$fc, xlab="", ylab="", main="", col=gray(0:100/100), axes=FALSE, zlim=zd11)
mtext("BEPT", side=2, line=0.3, cex=0.85, font=2)
image(outd11sift14$residue, xlab="", ylab="", main="", col=gray(0:100/100), axes=FALSE, zlim=zd11)




