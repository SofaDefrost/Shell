mt <- read.table("disp-tet.txt")
mc <- read.table("disp-comp.txt")

mt[,3] <- (mt[,3] - 0.015)*1000
mc[,3] <- (mc[,3] - 0.015)*1000


pdf("displacement.pdf", width=3, height=1.5, pointsize=8); # sizes in inches
p <- par(lwd=0.5, mar=c(4,4,0,0)+0.1);

i <- 1:200 # 5 seconds
plot(c(0,2),
     c(-0.285, 2.336), # make the range consistent
     #range(mt[,3], mc[,3]),
     type="n", xlab="Time [s]", ylab="Displacement [mm]")

lines(mt[i,1], mt[i,3], col=1)
lines(mc[i,1], mc[i,3], col=4)

dev.off()
par(p)
