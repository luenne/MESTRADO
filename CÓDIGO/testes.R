cluster = c('00053', '00086', '00339', '00719', '00996', 
'01052', '01189', '01238', '01264', '01347', '01478', 
'01831', '01836', '01877', '01933', '02035', '02104', 
'02137', '02182', '02186', '02249', '02298', '02301', 
'02433', '02440', '02447', '02469', '02490', '02752', 
'02759', '02789', '02827', '02899', '02907', '03031', 
'03112', '03176', '03229', '03459', '03565', '03631', 
'03691', '03742', '03898', '03907', '03915', '03975', 
'04023', '04048', '04100', '04376', '04404', '04405', 
'04409', '04458', '04470', '04479', '04537', '04641', 
'04672', '04681', '04703', '04710', '05039', '05206', 
'05325', '05359', '05447', '05535', '05717', '05780', 
'05859', '05908', '06070', '06173', '06175', '06184', 
'06207', '06233', '06256', '06261', '06264', '06286', 
'06392', '06447', '06475', '06506', '06508', '06547', 
'06723', '06758', '06841', '06861', '06924', '07204', 
'07217', '07395', '07435', '07520', '07584', '07703', 
'07775', '07815', '07837', '07975', '08022', '08049', 
'08173', '08219', '08291', '08338', '08710', '08720', 
'08721', '08738', '08742', '08975', '09061', '09075', 
'09132', '09144', '09148', '09153', '09157', '09162', 
'09176', '09177', '10001', '10004', '10006', '10008', 
'10009', '10010', '10013', '10014', '10015', '10016', 
'10017', '10018', '10019', '10020', '10021', '10022', 
'10023', '10024', '10025', '10026', '10027', '10028', 
'10029', '10030', '10031', '10032', '10033', '10034', 
'10035', '10036', '10037', '10038', '10039', '10040', 
'10041', '10042', '10043', '10044', '10045', '10046', 
'10047', '10048', '10049', '10050', '10051', '10052', 
'10053', '10054', '10055', '10056', '10058', '10059', 
'10060', '10062', '10063', '10064')

menor20 = 0
for(cl in 1:length(cluster)) {
	dir = paste("clusters_nosocs/cluster.gals.sel.shiftgap.iter.", cluster[cl], sep="")
	myInput = read.table(dir)

	inputT = data.frame(ra = myInput$V1, dec = myInput$V2, zspec = myInput$V8, Vel = myInput$V10,flag = myInput$V15)
	input <- inputT[which(inputT$flag == 1),names(inputT) %in% c("ra","dec","zspec","Vel","flag")]

	if(nrow(input) <= 20)
		menor20 = menor20 + 1
}

print(menor20)


# teste de normalidade

library(nortest)

cluster = c('00053', '00086', '00339', '00719', '00996', 
'01052', '01189', '01238', '01264', '01347', '01478', 
'01831', '01836', '01877', '01933', '02035', '02104', 
'02137', '02182', '02186', '02249', '02298', '02301', 
'02433', '02440', '02447', '02469', '02490', '02752', 
'02759', '02789', '02827', '02899', '02907', '03031', 
'03112', '03176', '03229', '03459', '03565', '03631', 
'03691', '03742', '03898', '03907', '03915', '03975', 
'04023', '04048', '04100', '04376', '04404', '04405', 
'04409', '04458', '04470', '04479', '04537', '04641', 
'04672', '04681', '04703', '04710', '05039', '05206', 
'05325', '05359', '05447', '05535', '05717', '05780', 
'05859', '05908', '06070', '06173', '06175', '06184', 
'06207', '06233', '06256', '06261', '06264', '06286', 
'06392', '06447', '06475', '06506', '06508', '06547', 
'06723', '06758', '06841', '06861', '06924', '07204', 
'07217', '07395', '07435', '07520', '07584', '07703', 
'07775', '07815', '07837', '07975', '08022', '08049', 
'08173', '08219', '08291', '08338', '08710', '08720', 
'08721', '08738', '08742', '08975', '09061', '09075', 
'09132', '09144', '09148', '09153', '09157', '09162', 
'09176', '09177', '10001', '10004', '10006', '10008', 
'10009', '10010', '10013', '10014', '10015', '10016', 
'10017', '10018', '10019', '10020', '10021', '10022', 
'10023', '10024', '10025', '10026', '10027', '10028', 
'10029', '10030', '10031', '10032', '10033', '10034', 
'10035', '10036', '10037', '10038', '10039', '10040', 
'10041', '10042', '10043', '10044', '10045', '10046', 
'10047', '10048', '10049', '10050', '10051', '10052', 
'10053', '10054', '10055', '10056', '10058', '10059', 
'10060', '10062', '10063', '10064')

ad_test = c()
for(cl in 1:length(cluster)) {
	dir = paste("clusters_nosocs/cluster.gals.sel.shiftgap.iter.", cluster[cl], sep="")
	myInput = read.table(dir)

	inputT = data.frame(Vel = myInput$V14,flag = myInput$V15, sigma = myInput$V17)
	input <- inputT[which(inputT$flag == 0),names(inputT) %in% c("Vel","flag", "sigma")]

	ad_test[cluster[cl]] = ad.test(input$Vel)$p.value
}
write.table(teste, file='arquivo.csv', sep=';', dec=',', row.names=FALSE)

clusters = which(ad_test < 0.1)



# definição de velocidade de rotação máxima eq. 24
myInput = read.table("info.sort")
input = data.frame(cl = myInput$V1, R200 = myInput$V27, M200 = myInput$V30)

G = 4.30e-9

input$M200 = input$M200 * 1e+14
omega = sqrt((G*input$M200)/(input$R200^3))


abline(h = 2.211725e-17,lty=5, lwd=2, col="red")

# teste de normalidade amostra SORIANO

library(nortest)

ad_test = c()
library(nortest)
myInput = read.table("rotacao/fort.8")
inputT = data.frame(ID = myInput$V1,zspec = myInput$V4, flag = myInput$V5)
for(ID in 1:max(inputT$ID)){
	input = inputT[which(inputT$ID == ID & inputT$flag == 1),names(inputT) %in% c("ID", "zspec")]
	velocity = input$zspec * 30000

	ad_test[ID] = ad.test(velocity)$p.value
}

write.table(ad_test, file='arquivo.csv', sep=';', dec=',', row.names=FALSE)



library(nortest)

ad_test = c()
library(nortest)
myInput = read.table("selec20.dat")
inputT = data.frame(ID = myInput$clusterID,Vel = myInput$Voff, flag = myInput$Flag)
for(ID in 1:max(inputT$ID)){
	input = inputT[which(inputT$ID == ID & inputT$flag == 0),names(inputT) %in% c("ID", "Vel")]

	ad_test[ID] = ad.test(input$Vel)$p.value
}

write.table(ad_test, file='arquivo.csv', sep=';', dec=',', row.names=FALSE)


# Corneta

dir = paste("clusters_nosocs/cluster.gals.sel.shiftgap.iter.10043", sep="")
myInput = read.table(dir)

inputT = data.frame(Voff = myInput$V14,flag = myInput$V15, radmpc = myInput$V13, zspec = myInput$V8)
cl01 <- inputT[which(inputT$flag == 0),names(inputT) %in% c("Voff", "radmpc", "zspec")]

rc=run_caustic(cl01$radmpc, cl01$Voff, median(cl01$zspec), H0 = 67.4)

 names(rc)

plot(cl01$radmpc,cl01$Voff,pch=20, ylab="Velocidade [km/s]", xlab="Distância ao centro [Mpc]")
plot(cl01$radmpc,cl01$Voff,pch=20, ylab="Velocidade [km/s]", xlab="Distância ao centro [Mpc]")
points(rc$x_range,rc$caustic_profile,type="l",lwd=2,col="darkorchid4")
points(rc$x_range,-1*rc$caustic_profile,type="l",lwd=2,col="darkorchid4")

# histograma

dir = paste("clusters_nosocs/cluster.gals.sel.shiftgap.iter.10043", sep="")
myInput = read.table(dir)

inputT = data.frame(Vel = myInput$V10,flag = myInput$V15)
input <- inputT[which(inputT$flag == 0),names(inputT) %in% c("Vel", "flag")]
input = input[order(input$Vel),]

h = hist(input$Vel, main = "Distribuição de Velocidades", 
 xlim = c(min(input$Vel),max(input$Vel)), breaks = 30)
plot(h$mids, h$counts, type="s",
xlab=expression(paste(Velocidades,"  ","(","km","/",s^{-1},")")), ylab="N")

xfit = seq(min(input$Vel), max(input$Vel), length=100)
yfit = dnorm(xfit,mean=mean(input$Vel), sd = sd(input$Vel))
yfit = yfit*diff(h$mids[1:2]*length(input$Vel))
lines(xfit,yfit,col="blue",lwd=2)
rug(input$Vel, ticksize=0.03, side=1, lwd=2, col="black")
