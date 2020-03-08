library(cosmoFns)

# aglomerados para leitura
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

# cluster = c('10062')
for(cl in 1:length(cluster)) {

	# leitura de arquivos para cada aglomerado
	dir = paste("clusters_nosocs/cluster.gals.sel.shiftgap.iter.", cluster[cl], sep="")
	myInput = read.table(dir)

	inputT = data.frame(ra = myInput$V1, dec = myInput$V2, zspec = myInput$V8, Vel = myInput$V10,flag = myInput$V15, r = myInput$V5)
	input <- inputT[which(inputT$flag == 0),names(inputT) %in% c("ra","dec","zspec","Vel","flag", "r")]

	if(nrow(input) > 20) {
		# cálculo da distância angular do sistema
		convrad = 0.0174533
		RA = input$ra*convrad
		DEC = input$dec*convrad

		decCl = median(DEC)
		raCl = median(RA)

		# Aglomerado no plano do céu
		mediumZ = median(input$zspec)

		groupdist = D.A(mediumZ, omega.m = 0.3, omega.lambda = 0.7, H.0 = 100)

		A = cos(DEC)*cos(RA - raCl)

		x = groupdist*(cos(DEC)*sin(RA - raCl))/(sin(decCl)*sin(DEC) + cos(decCl)*A)
		y = groupdist*(cos(decCl)*sin(DEC) - sin(decCl)*A)/(sin(decCl)*sin(DEC) + cos(decCl)*A)

		theta =  atan(y/x)

		vlos = input$Vel
		
		# velocidade e grau de rotação
		vRot = seq(-600, 600, length.out = nrow(input))
		vartheta = seq(0, pi, length.out = nrow(input))

		results_HL = data.frame()

		for(v in 1:length(vRot)) {
			vrot = vRot[v]
			chiq2 = vector()
			for(t in 1:length(vartheta)){
				theta0 = vartheta[t]
				chiq2[t] =  sum((vrot^2  *  (sin(theta - theta0))^2)/ (vlos + vrot * sin(theta - theta0)))
			}

			chi2=min(chiq2)
			min_theta0 = vartheta[which.min(chiq2)] ## rad
			min_theta0_deg = min_theta0 * 57.29577 ## degrees

			resvec=cbind(vrot,chi2,min_theta0,min_theta0_deg)
			results_HL = rbind(results_HL,resvec)
		}

		
		res_final = results_HL[which.min(results_HL$chi2),]
		vrot_final =  res_final$vrot
		theta0_final = res_final$min_theta0
		theta0_final_deg = theta0_final*57.29577
		print(cluster[cl])
		print(res_final)
		# write.table(paste(res_final), file='arquivo3.csv', sep=';', dec=',', row.names=FALSE, append=TRUE)
	}
}