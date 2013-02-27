a1 <- data.matrix(read.csv("aspiration-simple.csv"))
a2 <- data.matrix(read.csv("aspiration-composite.csv"))

a1 <- a1[1:21,c(4,3)]; a1[,2] <- a1[,2] - 0.03
a2 <- a2[1:21,c(4,3)]; a2[,2] <- a2[,2] - 0.03

pdf("aspiration.pdf", width=4, height=2, pointsize=9); # sizes in inches
p <- par(lwd=2, mar=c(4,4,1,1)+0.1);

plot(range(c(a1[,1], a2[,1])), range(c(a1[,2], a2[,2])), #ylim=c(0, 0.004), 
     type="n", xlab="x-coordinate", ylab="Displacement")

lines(a1, col="black")
lines(a2, col="blue")

abline(h=0, lty="dotted");

#legend("topright", c("A","B", "C", "D"), pch=c(1,4,5,6), col="black");

dev.off()
par(p);
