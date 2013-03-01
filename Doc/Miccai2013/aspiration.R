a1 <- data.matrix(read.csv("aspiration-simple.csv"))
a2 <- data.matrix(read.csv("aspiration-composite.csv"))

a1 <- a1[1:31,c(4,3)]; a1[,2] <- a1[,2] - 0.015; a1[,2] <- a1[,2] * 1000
a2 <- a2[1:33,c(4,3)]; a2[,2] <- a2[,2] - 0.015; a2[,2] <- a2[,2] * 1000

pdf("aspiration.pdf", width=3, height=1.5, pointsize=8); # sizes in inches
p <- par(lwd=0.5, mar=c(4,4,0,0)+0.1);

plot(range(c(a1[,1], a2[,1])), range(c(a1[,2], a2[,2])), #ylim=c(0, 0.004), 
     type="n", xlab="x-coordinate", ylab="Displacement [mm]")

lines(a1, col="black")
lines(a2, col="blue")

abline(h=0, lty="dotted");

#legend("topright", c("A","B", "C", "D"), pch=c(1,4,5,6), col="black");

dev.off()
par(p);
