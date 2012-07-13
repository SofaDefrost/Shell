#
# Graph: Frames per seconds with and withouth frame fixing
#

fps <- read.csv("fps.csv");

pdf("fps.pdf", width=7, height=4);
p <- par(lwd=2, mar=c(4,4,1,1)+0.1);

plot(range(fps$elements), range(c(fps$fpsBase, fps$fpsPolar, fps$fpsMap1, fps$fpsMap2)),
    type="n", xlab="Elements", ylab="FPS", log="y");

lines(fps$elements, fps$fpsBase, type="b", col="black", pch=1);
lines(fps$elements, fps$fpsPolar, type="b", col="black", pch=4);
lines(fps$elements, fps$fpsMap2, type="b", col="black", pch=6);
lines(fps$elements, fps$fpsMap1, type="b", col="black", pch=5);

abline(h=25, lty="dotted");

legend("topright", c("A","B", "C", "D"), pch=c(1,4,5,6), col="black");

dev.off()
par(p);
