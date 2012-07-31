#
# Graph: modified Scordelis-Lo roof
#

n2e <- function(x) (2*4*x*x); # size to element count

hem <- read.csv("hemisphere.csv");

# Compute radial displacement
d <- apply(cbind(hem$x,hem$y,hem$z), 1, function(x) (sqrt(sum(x*x)))) - 10
hem <- data.frame(hem, d=d);

# Plot to PDF
pdf("hemisphere.pdf", width=7, height=4);
p <- par(lwd=2, mar=c(4,4,1,1)+0.1);

plot(c(1,3200), c(0,0.2), type="n", xlab="Elements", ylab="Displacement");
abline(h=axTicks(2), lty="dotted", lwd=1); # Grid lines

n <- "optdkt";
lines(n2e(hem$size[hem$name == n]), hem$d[hem$name == n], type="b", col="black", pch=1);

n<- "bezier2";
lines(n2e(hem$size[hem$name == n]), hem$d[hem$name == n], type="b", col="black", pch=4);

n<- "triangular";
lines(n2e(hem$size[hem$name == n]), hem$d[hem$name == n], type="b", col="black", pch=5);

legend("bottomright", c("DKT","BSH","PSH"), pch=c(1,4,5), col="black", bg="white");

dev.off()
par(p);
