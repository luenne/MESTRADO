library(astro)
library(ellipse)
library(cosmoFns)
library(cramer)
library(Hotelling)

 dist_2_points <- function(x1, x2) {
     return(sqrt(sum((x1 - x2)^2)))    
 }

  # função distância entre um ponto e uma reta
 dist_point_line <- function(a, slope, intercept) {
    b = c(1, intercept+slope)
    c = c(-intercept/slope,0)      
    v1 <- b - c
    v2 <- a - b
    m <- cbind(v1,v2)
    
    return(det(m))/sqrt(sum(v1*v1))
}

# função omega
omegaV <- function(a_point, b_point, center, R200, M200) {
	G = 4.30e-9 # mpc Msun km/s

	a = dist_2_points(a_point, center)
	b = dist_2_points(b_point, center)

	e = 1 - (b/a)
	velocity = sqrt((1/R200^4) * G * M200 * a * (1 - e^2))

	return(velocity/(R200*3.08e19))
}

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


sample = c()
rot = matrix(nrow = 1000, ncol = 183)

for(cl in 1:length(cluster)) {
	dir = paste("clusters_nosocs/cluster.gals.sel.shiftgap.iter.", cluster[cl], sep="")
	myInput = read.table(dir)
	
	inputT = data.frame(ra = myInput$V1, dec = myInput$V2, zspec = myInput$V8, Vel = myInput$V10,flag = myInput$V15, r = myInput$V5)

	input <- inputT[which(inputT$flag == 0),names(inputT) %in% c("ra","dec","zspec","Vel","flag", "r")]

	for(it in 1:1000) {
		if(nrow(input) <= 20)
			sample[it] = "MENOR"
		else {
			input = input[order(input$Vel),]

			g = c()
			w = c()
			weiGap = c()
			N = nrow(input) # número de objetos do aglomerado
			lowerS = round(N/4)
			upS = round(3*N/4)

			for(i in 1:(N-1)) {
				g[i] = input$Vel[i+1] - input$Vel[i] #gap
				w[i] = i*(N-i) #peso
				weiGap[i] = sqrt(w[i]*g[i])
			}

			MM = (2/N)*sum(weiGap[lowerS:upS])
			weiGap = weiGap/MM

			# Considera apenas gaps > 2.25
			gap_positions = which(weiGap > 2.25)

			# gap não encontrado, calcula mediana
			if(length(gap_positions) == 0){
				medVel = input$Vel[which(input$Vel >= median(input$Vel))]
				gap_positions = medVel[1]
			}

			pos_gap = which.max(weiGap)

			input['Flag'] = 0

			# Valores à esquerda do maior gap recebe em Flag o valor -1 e à direita 1
			for(i in 1:N) {
				if(pos_gap >= i) 
					input$Flag[i] = -1
				else
					input$Flag[i] = 1
			}

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

			# Insere os valores de x e y no dataframe
			input['Vx'] = x
			input['Vy'] = y

			eli = ellipse(cor(x,y),scale=c(sd(x),sd(y)),centre=c(mean(x),mean(y)),level=0.95, npoints=round(length(x)*0.95))

			# Calcula o centro da elipse
			center_elip = c(mean(eli[,1]), mean(eli[,2]))

			distance = numeric(0)
			for(i in 1:nrow(eli)) {
				distance[i] = dist_2_points(center_elip, eli[i,])
			}

			# Distância máxima ao centro da elipse (vertical)
			a = distance[which.max(distance)]
			a_point = eli[which.max(distance),]

			# Distância máxima ao centro da elipse (horizontal)
			b = distance[which.min(distance)]
			b_point = eli[which.min(distance),]

			xx = c(center_elip[1], a_point[1])
			yy = c(center_elip[2],a_point[2])

			fitline = lm(yy ~ xx)
			alpha = fitline$coef[2]	# inclinação
			beta = fitline$coef[1]	# intersecção

			dx = seq(-2,2,0.02)
			ypred = dx*alpha + beta

			pnts = cbind(x,y)
			distAxis = apply(pnts, 1, function(x) dist_point_line(as.numeric(x[1:2]), slope = alpha, intercept = beta))

			# Insere as distâncias na última coluna do dataframe
			input['DistA'] = distAxis
		 	
		 	rotacao = FALSE # considerando que não haja uma possível rotação

			# Todos os pontos do gráfico
			a1 = subset(input[c("Vx","Vy")], input$Flag == -1)
			a2 = subset(input[c("Vx","Vy")], input$Flag == 1)
			a1 = matrix(c(a1$Vx, a1$Vy), ncol = 2)
			a2 = matrix(c(a2$Vx, a2$Vy), ncol = 2)

			if(length(a1) > 0 && length(a2) > 0) {
				cramer = cramer.test(a1,a2)
				
				hotelling = hotelling.test(a1,a2)
				
				if((!is.na(cramer$p.value) && cramer$p.value < 0.05) || (!is.na(hotelling$pval) && hotelling$pval < 0.05))
					rotacao = TRUE
			}
			# Todos os pontos acima do eixo principal
			a1 = subset(input[c("Vx","Vy")], input$Flag == -1 & input$DistA >= 0)
			a2 = subset(input[c("Vx","Vy")], input$Flag == 1 & input$DistA >= 0)
			a1 = matrix(c(a1$Vx,a1$Vy), ncol=2)
			a2 = matrix(c(a2$Vx,a2$Vy), ncol=2)

			if(length(a1) > 0 && length(a2) > 0) {
				cramer = cramer.test(a1,a2)
				
				hotelling = hotelling.test(a1,a2)
				
				if((!is.na(cramer$p.value) && cramer$p.value < 0.05) || (!is.na(hotelling$pval) && hotelling$pval < 0.05))
					rotacao = TRUE
			}

			# Todos os pontos abaixo do eixo principal
			a1 = subset(input[c("Vx","Vy")], input$Flag == -1 & input$DistA < 0)
			a2 = subset(input[c("Vx","Vy")], input$Flag == 1 & input$DistA < 0)
			a1 = matrix(c(a1$Vx,a1$Vy), ncol=2)
			a2 = matrix(c(a2$Vx,a2$Vy), ncol=2)

			if(length(a1) > 0 && length(a2) > 0) {
				cramer = cramer.test(a1,a2)
				
				hotelling = hotelling.test(a1,a2)

				if((!is.na(cramer$p.value) && cramer$p.value < 0.05) || (!is.na(hotelling$pval) && hotelling$pval < 0.05))
					rotacao = TRUE
			}

			sample[it] = rotacao 

		} # final > 20

	} # final 1000 simulações

	print(cluster[cl])
	rot[ , cl] = sample
	
} # final aglomerados
