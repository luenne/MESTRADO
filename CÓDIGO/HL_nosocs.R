library(cosmoFns)

cluster = c('00339')

for(cl in 1:length(cluster)) {

	# Leitura de arquivo
	dir = paste("clusters_nosocs/cluster.gals.sel.shiftgap.iter.", cluster[cl], sep="")
	myInput = read.table(dir)
	
	inputT = data.frame(ra = myInput$V1, dec = myInput$V2, vlos = myInput$V14, zspec = myInput$V8, flag = myInput$V15)

	input <- inputT[which(inputT$flag == 0),names(inputT) %in% c("ra","dec", "vlos","zspec","flag")]

	# --- CALCULA A DISTÂNCIA ANGULAR DO SISTEMA ---

	# Conversão para radiano
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

	theta = atan(y/x)

	vRot = seq(-1184.30, 1085.40, length.out = nrow(input))
	vartheta = seq(0, pi, length.out = nrow(input))
	
	results_HL = data.frame()

	for(v in 1:length(vRot)) {
		vrot = vRot[v]
		chiq2 = vector()
		for(t in 1:length(vartheta)){
			theta0 = vartheta[t]
			chiq2[t] =  sum((vrot^2  *  (sin(theta - theta0))^2)/ (input$vlos + vrot * sin(theta - theta0)))
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
	
}
