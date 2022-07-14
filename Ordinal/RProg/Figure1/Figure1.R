#####################################################################################################################
#####################################################################################################################
#####################################################################################################################
#####################################################################################################################
#####################################################################################################################
#####################################################################################################################
#tiff("C:\\Bayes_PXDA_MVP_Ordinal\\CRPrograms-SimuStudy\\PX-Analyses\\SimuStudyResults\\Figure_1.tif", width = 4800, height = 4800,res = 800, compression = "lzw")
par(mfrow=c(4,4), mar=c(4,4,1,0.1))


###par(mfrow = c(4, 4))
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################


plot(pb150ID, type = "l", frame = FALSE, pch = 19, col = "black", xlab = bquote(beta[1]), ylab = "n=50 PID", lty = 1, lwd = 3)
lines(sb150ID, pch = 18, col = "blue", type = "l", lty = 5, lwd = 3)
lines(wb150ID, pch = 18, col = "purple", type = "l", lty = 4, lwd = 3)


plot(pg1150ID, type = "l", frame = FALSE, pch = 19, col = "black", xlab = bquote(gamma[11]), ylab = " ", lty = 1, lwd = 3)
lines(sg1150ID, pch = 18, col = "blue", type = "l", lty = 5, lwd = 3)
lines(wg1150ID, pch = 18, col = "purple", type = "l", lty = 4, lwd = 3)


plot(pr2350ID, type = "l", frame = FALSE, pch = 19, col = "black", xlab = bquote(r[23]), ylab = " ", lty = 1, lwd = 3)
lines(sr2350ID, pch = 18, col = "blue", type = "l", lty = 5, lwd = 3)
lines(wr2350ID, pch = 18, col = "purple", type = "l", lty = 4, lwd = 3)


plot(pr1550ID, type = "l", frame = FALSE, pch = 19, col = "black", xlab = bquote(r[15]), ylab = " ", lty = 1, lwd = 3)
lines(sr1550ID, pch = 18, col = "blue", type = "l", lty = 5, lwd = 3)
lines(wr1550ID, pch = 18, col = "purple", type = "l", lty = 4, lwd = 3)

#######################################################################################################################


plot(pb150CS, type = "l", frame = FALSE, pch = 19, col = "black", xlab = bquote(beta[1]), ylab = "n=50 PCS", lty = 1, lwd = 3)
lines(sb150CS, pch = 18, col = "blue", type = "l", lty = 5, lwd = 3)
lines(wb150CS, pch = 18, col = "purple", type = "l", lty = 4, lwd = 3)


plot(pg1150CS, type = "l", frame = FALSE, pch = 19, col = "black", xlab = bquote(gamma[11]), ylab = " ", lty = 1, lwd = 3)
lines(sg1150CS, pch = 18, col = "blue", type = "l", lty = 5, lwd = 3)
lines(wg1150CS, pch = 18, col = "purple", type = "l", lty = 4, lwd = 3)


plot(pr2350CS, type = "l", frame = FALSE, pch = 19, col = "black", xlab = bquote(r[23]), ylab = " ", lty = 1, lwd = 3)
lines(sr2350CS, pch = 18, col = "blue", type = "l", lty = 5, lwd = 3)
lines(wr2350CS, pch = 18, col = "purple", type = "l", lty = 4, lwd = 3)


plot(pr1550CS, type = "l", frame = FALSE, pch = 19, col = "black", xlab = bquote(r[15]), ylab = " ", lty = 1, lwd = 3)
lines(sr1550CS, pch = 18, col = "blue", type = "l", lty = 5, lwd = 3)
lines(wr1550CS, pch = 18, col = "purple", type = "l", lty = 4, lwd = 3)

#######################################################################################################################



plot(pb1500ID, type = "l", frame = FALSE, pch = 19, col = "black", xlab = bquote(beta[1]), ylab = "n=500 PID", lty = 1, lwd = 3)
lines(sb1500ID, pch = 18, col = "blue", type = "l", lty = 5, lwd = 3)
lines(wb1500ID, pch = 18, col = "purple", type = "l", lty = 4, lwd = 3)


plot(pg11500ID, type = "l", frame = FALSE, pch = 19, col = "black", xlab = bquote(gamma[11]), ylab = " ", lty = 1, lwd = 3)
lines(sg11500ID, pch = 18, col = "blue", type = "l", lty = 5, lwd = 3)
lines(wg11500ID, pch = 18, col = "purple", type = "l", lty = 4, lwd = 3)


plot(pr23500ID, type = "l", frame = FALSE, pch = 19, col = "black", xlab = bquote(r[23]), ylab = " ", lty = 1, lwd = 3)
lines(sr23500ID, pch = 18, col = "blue", type = "l", lty = 5, lwd = 3)
lines(wr23500ID, pch = 18, col = "purple", type = "l", lty = 4, lwd = 3)


plot(pr15500ID, type = "l", frame = FALSE, pch = 19, col = "black", xlab = bquote(r[15]), ylab = " ", lty = 1, lwd = 3)
lines(sr15500ID, pch = 18, col = "blue", type = "l", lty = 5, lwd = 3)
lines(wr15500ID, pch = 18, col = "purple", type = "l", lty = 4, lwd = 3)


#######################################################################################################################



plot(pb1500CS, type = "l", frame = FALSE, pch = 19, col = "black", xlab = bquote(beta[1]), ylab = "n=500 PCS", lty = 1, lwd = 3)
lines(sb1500CS, pch = 18, col = "blue", type = "l", lty = 5, lwd = 3)
lines(wb1500CS, pch = 18, col = "purple", type = "l", lty = 4, lwd = 3)


plot(pg11500CS, type = "l", frame = FALSE, pch = 19, col = "black", xlab = bquote(gamma[11]), ylab = " ", lty = 1, lwd = 3)
lines(sg11500CS, pch = 18, col = "blue", type = "l", lty = 5, lwd = 3)
lines(wg11500CS, pch = 18, col = "purple", type = "l", lty = 4, lwd = 3)


plot(pr23500CS, type = "l", frame = FALSE, pch = 19, col = "black", xlab = bquote(r[23]), ylab = " ", lty = 1, lwd = 3)
lines(sr23500CS, pch = 18, col = "blue", type = "l", lty = 5, lwd = 3)
lines(wr23500CS, pch = 18, col = "purple", type = "l", lty = 4, lwd = 3)


plot(pr15500CS, type = "l", frame = FALSE, pch = 19, col = "black", xlab = bquote(r[15]), ylab = " ", lty = 1, lwd = 3)
lines(sr15500CS, pch = 18, col = "blue", type = "l", lty = 5, lwd = 3)
lines(wr15500CS, pch = 18, col = "purple", type = "l", lty = 4, lwd = 3)



#######################################################################################################################
###########################################################################################################################
###########################################################################################################################
###########################################################################################################################
#dev.off() 
###########################################################################################################################
###########################################################################################################################
###########################################################################################################################
###########################################################################################################################


############################################################################################################################
############################################################################################################################
############################################################################################################################
############################################################################################################################












