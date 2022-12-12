# plot 1: Proportion

proportion10_30 <- read.csv("./Dropbox/Research_Documents/AP/Collaboration/R-CDDU/julia_code/proportion_10_30.csv", header = FALSE)
png(file = "./Dropbox/Research_Documents/AP/Collaboration/R-CDDU/julia_code/proportion_10_30.png", width= 10,height=6,units = 'in',res = 300);
par(mar = c(5,5,2,2));
plot(proportion10_30$V1*0.01, proportion10_30$V2*100,type = "l",ylim=range(c(0,60)),xlab = expression(paste("Uncertainty Budget ", Gamma)), ylab = "Proportion of Commitment (%)", 
     main = "Commitment Proportion of DR Resource Types", col = "#377EB8", lwd = 3, lty = 2, cex.main = 2, cex.lab = 2, cex.axis = 2);
lines(proportion10_30$V1*0.01, proportion10_30$V3*100, col = "#377EB8", lwd = 3, lty = 1)
lines(proportion10_30$V1*0.01, proportion10_30$V4*100, col = "#377EB8", lwd = 3, lty = 4)
legend("bottomright",c("Type-A","Type-B","Type-C"), col = c("#377EB8","#377EB8","#377EB8"), lty = c(2,1,4), lwd = 3, cex = 1.3)
dev.off();

proportion10_80 <- read.csv("./Dropbox/Research_Documents/AP/Collaboration/R-CDDU/julia_code/proportion_10_80.csv", header = FALSE)
png(file = "./Dropbox/Research_Documents/AP/Collaboration/R-CDDU/julia_code/proportion_10_80.png", width= 10,height=6,units = 'in',res = 300);
par(mar = c(5,5,2,2));
plot(proportion10_80$V1*0.01, proportion10_80$V2*100,type = "l",ylim=range(c(0,60)),xlab = expression(paste("Uncertainty Budget ", Gamma)), ylab = "Proportion of Commitment (%)", 
     main = "Commitment Proportion of DR Resource Types", col = "#377EB8", lwd = 3, lty = 2, cex.main = 2, cex.lab = 2, cex.axis = 2);
lines(proportion10_80$V1*0.01, proportion10_80$V3*100, col = "#377EB8", lwd = 3, lty = 1)
lines(proportion10_80$V1*0.01, proportion10_80$V4*100, col = "#377EB8", lwd = 3, lty = 4)
legend("bottomright",c("Type-A","Type-B","Type-C"), col = c("#377EB8","#377EB8","#377EB8"), lty = c(2,1,4), lwd = 3, cex = 1.3)
dev.off();

proportion30_30 <- read.csv("./Dropbox/Research_Documents/AP/Collaboration/R-CDDU/julia_code/proportion_30_30.csv", header = FALSE)
png(file = "./Dropbox/Research_Documents/AP/Collaboration/R-CDDU/julia_code/proportion_30_30.png", width= 10,height=6,units = 'in',res = 300);
par(mar = c(5,5,2,2));
plot(proportion30_30$V1*0.01, proportion30_30$V2*100,type = "l",ylim=range(c(0,60)),xlab = expression(paste("Uncertainty Budget ", Gamma)), ylab = "Proportion of Commitment (%)", 
     main = "Commitment Proportion of DR Resource Types", col = "#377EB8", lwd = 3, lty = 2, cex.main = 2, cex.lab = 2, cex.axis = 2);
lines(proportion30_30$V1*0.01, proportion30_30$V3*100, col = "#377EB8", lwd = 3, lty = 1)
lines(proportion30_30$V1*0.01, proportion30_30$V4*100, col = "#377EB8", lwd = 3, lty = 4)
legend("bottomright",c("Type-A","Type-B","Type-C"), col = c("#377EB8","#377EB8","#377EB8"), lty = c(2,1,4), lwd = 3, cex = 1.3)
dev.off();

proportion30_80 <- read.csv("./Dropbox/Research_Documents/AP/Collaboration/R-CDDU/julia_code/proportion_30_80.csv", header = FALSE)
png(file = "./Dropbox/Research_Documents/AP/Collaboration/R-CDDU/julia_code/proportion_30_80.png", width= 10,height=6,units = 'in',res = 300);
par(mar = c(5,5,2,2));
plot(proportion30_80$V1*0.01, proportion30_80$V2*100,type = "l",ylim=range(c(0,60)),xlab = expression(paste("Uncertainty Budget ", Gamma)), ylab = "Proportion of Commitment (%)", 
     main = "Commitment Proportion of DR Resource Types", col = "#377EB8", lwd = 3, lty = 2, cex.main = 2, cex.lab = 2, cex.axis = 2);
lines(proportion30_80$V1*0.01, proportion30_80$V3*100, col = "#377EB8", lwd = 3, lty = 1)
lines(proportion30_80$V1*0.01, proportion30_80$V4*100, col = "#377EB8", lwd = 3, lty = 4)
legend("bottomright",c("Type-A","Type-B","Type-C"), col = c("#377EB8","#377EB8","#377EB8"), lty = c(2,1,4), lwd = 3, cex = 1.3)
dev.off();


#=============================================================================================
# plot 2: 
DRamount10_30 <- read.csv("./Dropbox/Research_Documents/AP/Collaboration/R-CDDU/julia_code/DRamount_10_30.csv", header = FALSE)
png(file = "./Dropbox/Research_Documents/AP/Collaboration/R-CDDU/julia_code/DRamount_10_30.png", width= 10,height=6,units = 'in',res = 300);
par(mar = c(5,5,2,2));
plot(DRamount10_30$V1, DRamount10_30$V2,type = "l",ylim=range(c(0,7000)), xlab = "Time", ylab = "Commitment Amount", 
     main = "Total Commitment Amount vs. Time", col = "#377EB8", lwd = 3, lty = 2, cex.main = 2, cex.lab = 2, cex.axis = 2);
lines(DRamount10_30$V1, DRamount10_30$V3, col = "#377EB8", lwd = 3, lty = 1)
lines(DRamount10_30$V1, DRamount10_30$V4, col = "#E41A1C", lwd = 3, lty = 1)
lines(DRamount10_30$V1, DRamount10_30$V5, col = "#E41A1C", lwd = 3, lty = 2)
lines(DRamount10_30$V1, DRamount10_30$V6, col = "#E41A1C", lwd = 3, lty = 4)

legend("topleft",c("Nominal Demand","Deterministic",expression(paste("Robust Sol ",Gamma, "= 0.01")),expression(paste("Robust Sol ",Gamma, "= 0.05")),expression(paste("Robust Sol ",Gamma, "= 0.1"))), col = c("#377EB8","#377EB8","#E41A1C","#E41A1C","#E41A1C"), lty = c(2,1,1,2,4), lwd = 3, cex = 1.3)
dev.off();

DRamount10_80 <- read.csv("./Dropbox/Research_Documents/AP/Collaboration/R-CDDU/julia_code/DRamount_10_80.csv", header = FALSE)
png(file = "./Dropbox/Research_Documents/AP/Collaboration/R-CDDU/julia_code/DRamount_10_80.png", width= 10,height=6,units = 'in',res = 300);
par(mar = c(5,5,2,2));
plot(DRamount10_80$V1, DRamount10_80$V2,type = "l",ylim=range(c(0,16000)), xlab = "Time", ylab = "Commitment Amount", 
     main = "Total Commitment Amount vs. Time", col = "#377EB8", lwd = 3, lty = 2, cex.main = 2, cex.lab = 2, cex.axis = 2);
lines(DRamount10_80$V1, DRamount10_80$V3, col = "#377EB8", lwd = 3, lty = 1)
lines(DRamount10_80$V1, DRamount10_80$V4, col = "#E41A1C", lwd = 3, lty = 1)
lines(DRamount10_80$V1, DRamount10_80$V5, col = "#E41A1C", lwd = 3, lty = 2)
lines(DRamount10_80$V1, DRamount10_80$V6, col = "#E41A1C", lwd = 3, lty = 4)

legend("topleft",c("Nominal Demand","Deterministic",expression(paste("Robust Sol ",Gamma, "= 0.01")),expression(paste("Robust Sol ",Gamma, "= 0.05")),expression(paste("Robust Sol ",Gamma, "= 0.1"))), col = c("#377EB8","#377EB8","#E41A1C","#E41A1C","#E41A1C"), lty = c(2,1,1,2,4), lwd = 3, cex = 1.3)
dev.off();

DRamount30_30 <- read.csv("./Dropbox/Research_Documents/AP/Collaboration/R-CDDU/julia_code/DRamount_30_30.csv", header = FALSE)
png(file = "./Dropbox/Research_Documents/AP/Collaboration/R-CDDU/julia_code/DRamount_30_30.png", width= 10,height=6,units = 'in',res = 300);
par(mar = c(5,5,2,2));
plot(DRamount30_30$V1, DRamount30_30$V2,type = "l",ylim=range(c(0,7000)), xlab = "Time", ylab = "Commitment Amount", 
     main = "Total Commitment Amount vs. Time", col = "#377EB8", lwd = 3, lty = 2, cex.main = 2, cex.lab = 2, cex.axis = 2);
lines(DRamount30_30$V1, DRamount30_30$V3, col = "#377EB8", lwd = 3, lty = 1)
lines(DRamount30_30$V1, DRamount30_30$V4, col = "#E41A1C", lwd = 3, lty = 1)
lines(DRamount30_30$V1, DRamount30_30$V5, col = "#E41A1C", lwd = 3, lty = 2)
lines(DRamount30_30$V1, DRamount30_30$V6, col = "#E41A1C", lwd = 3, lty = 4)

legend("topleft",c("Nominal Demand","Deterministic",expression(paste("Robust Sol ",Gamma, "= 0.01")),expression(paste("Robust Sol ",Gamma, "= 0.05")),expression(paste("Robust Sol ",Gamma, "= 0.1"))), col = c("#377EB8","#377EB8","#E41A1C","#E41A1C","#E41A1C"), lty = c(2,1,1,2,4), lwd = 3, cex = 1.3)
dev.off();

DRamount30_80 <- read.csv("./Dropbox/Research_Documents/AP/Collaboration/R-CDDU/julia_code/DRamount_30_80.csv", header = FALSE)
png(file = "./Dropbox/Research_Documents/AP/Collaboration/R-CDDU/julia_code/DRamount_30_80.png", width= 10,height=6,units = 'in',res = 300);
par(mar = c(5,5,2,2));
plot(DRamount30_80$V1, DRamount30_80$V2,type = "l",ylim=range(c(0,16000)), xlab = "Time", ylab = "Commitment Amount", 
     main = "Total Commitment Amount vs. Time", col = "#377EB8", lwd = 3, lty = 2, cex.main = 2, cex.lab = 2, cex.axis = 2);
lines(DRamount30_80$V1, DRamount30_80$V3, col = "#377EB8", lwd = 3, lty = 1)
lines(DRamount30_80$V1, DRamount30_80$V4, col = "#E41A1C", lwd = 3, lty = 1)
lines(DRamount30_80$V1, DRamount30_80$V5, col = "#E41A1C", lwd = 3, lty = 2)
lines(DRamount30_80$V1, DRamount30_80$V6, col = "#E41A1C", lwd = 3, lty = 4)

legend("topleft",c("Nominal Demand","Deterministic",expression(paste("Robust Sol ",Gamma, "= 0.01")),expression(paste("Robust Sol ",Gamma, "= 0.05")),expression(paste("Robust Sol ",Gamma, "= 0.1"))), col = c("#377EB8","#377EB8","#E41A1C","#E41A1C","#E41A1C"), lty = c(2,1,1,2,4), lwd = 3, cex = 1.3)
dev.off();

#=============================================================================================
# plot 3:
simu10_30 <- read.csv("./Dropbox/Research_Documents/AP/Collaboration/R-CDDU/julia_code/simu_10_30.csv", header = FALSE)
png(file = "./Dropbox/Research_Documents/AP/Collaboration/R-CDDU/julia_code/simu_10_30.png", width= 10, height=6, units = 'in',res = 300);
par(mar = c(5,5,2,2));
plot(simu10_30$V1*0.01, simu10_30$V2/100000,ylim=range(c(-0.5,5)),xlab = expression(paste("Uncertainty Budget ", Gamma)), ylab = "",
     main = expression(paste("Simulated Profit vs. ",Gamma)), pch = 19, col = "#377EB8", cex.main = 2, cex.lab = 2, cex.axis = 2);
lines(simu10_30$V1*0.01, simu10_30$V2/100000, col = "#377EB8", lwd = 3, lty = 1)
title( ylab = expression(paste("Simulated Profit ($10"^"5",")")), line = 2.2, cex.lab = 2)
dev.off();

simu10_80 <- read.csv("./Dropbox/Research_Documents/AP/Collaboration/R-CDDU/julia_code/simu_10_80.csv", header = FALSE)
png(file = "./Dropbox/Research_Documents/AP/Collaboration/R-CDDU/julia_code/simu_10_80.png", width= 10, height=6, units = 'in',res = 300);
par(mar = c(5,5,2,2));
plot(simu10_80$V1*0.01, simu10_80$V2/100000,ylim=range(c(-0.5,5)),xlab = expression(paste("Uncertainty Budget ", Gamma)), ylab = "", 
     main = expression(paste("Simulated Profit vs. ",Gamma)), pch = 19, col = "#377EB8", cex.main = 2, cex.lab = 2, cex.axis = 2);
lines(simu10_80$V1*0.01, simu10_80$V2/100000, col = "#377EB8", lwd = 3, lty = 1)
title( ylab = expression(paste("Simulated Profit ($10"^"5",")")), line = 2.2, cex.lab = 2)
dev.off();

simu30_30 <- read.csv("./Dropbox/Research_Documents/AP/Collaboration/R-CDDU/julia_code/simu_30_30.csv", header = FALSE)
png(file = "./Dropbox/Research_Documents/AP/Collaboration/R-CDDU/julia_code/simu_30_30.png", width= 10, height=6, units = 'in',res = 300);
par(mar = c(5,5,2,2));
plot(simu30_30$V1*0.01, simu30_30$V2/100000,ylim=range(c(-0.5,5)),xlab = expression(paste("Uncertainty Budget ", Gamma)), ylab = "", 
     main = expression(paste("Simulated Profit vs. ",Gamma)), pch = 19, col = "#377EB8", cex.main = 2, cex.lab = 2, cex.axis = 2);
lines(simu30_30$V1*0.01, simu30_30$V2/100000, col = "#377EB8", lwd = 3, lty = 1)
title( ylab = expression(paste("Simulated Profit ($10"^"5",")")), line = 2.2, cex.lab = 2)
dev.off();

simu30_80 <- read.csv("./Dropbox/Research_Documents/AP/Collaboration/R-CDDU/julia_code/simu_30_80.csv", header = FALSE)
png(file = "./Dropbox/Research_Documents/AP/Collaboration/R-CDDU/julia_code/simu_30_80.png", width= 10, height=6, units = 'in',res = 300);
par(mar = c(5,5,2,2));
plot(simu30_80$V1*0.01, simu30_80$V2/100000,ylim=range(c(-0.5,5)),xlab = expression(paste("Uncertainty Budget ", Gamma)), ylab = "", 
     main = expression(paste("Simulated Profit vs. ",Gamma)), pch = 19, col = "#377EB8", cex.main = 2, cex.lab = 2, cex.axis = 2);
lines(simu30_80$V1*0.01, simu30_80$V2/100000, col = "#377EB8", lwd = 3, lty = 1)
title( ylab = expression(paste("Simulated Profit ($10"^"5",")")), line = 2.2, cex.lab = 2)
dev.off();

