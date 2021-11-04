#
# Graph: modified Scordelis-Lo roof
#

n2e <- function(x) (2*x*x); # size to element count

lo <- read.csv("roof.csv");

# Compute displacement
lo$z <- 19.1511 - lo$z;

# Plot to PDF
pdf("roof.pdf", width=7, height=4);
p <- par(lwd=3, mar=c(4,4,1,1)+0.1);

plot(c(1,3200), c(0,0.5), type="n", xlab="Elements", ylab="Displacement");
abline(h=axTicks(2), lty="dotted", lwd=1); # Grid lines

n <- "optdkt";
lines(n2e(lo$size[lo$name == n]), lo$z[lo$name == n], type="b", col="black", pch=1);

n<- "bezier2";
lines(n2e(lo$size[lo$name == n]), lo$z[lo$name == n], type="b", col="black", pch=4);

n<- "triangular";
lines(n2e(lo$size[lo$name == n]), lo$z[lo$name == n], type="b", col="black", pch=5);

legend("topright", c("DKT","BSH", "PSH"), pch=c(1,4,5), col=c("black","black","black"), bg="white");

dev.off()
par(p);
