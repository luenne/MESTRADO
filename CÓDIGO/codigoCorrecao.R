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

# cluster = c('00053', '00086', '00339', '00719', '00996', 
# '01052', '01189', '01238', '01264', '01347', '01478', 
# '01831', '01836', '01877', '01933', '02035', '02104', 
# '02137', '02182', '02186', '02249', '02298', '02301', 
# '02433', '02440', '02447', '02469', '02490', '02752', 
# '02759', '02789', '02827', '02899', '02907', '03031', 
# '03112', '03176', '03229', '03459', '03565', '03631', 
# '03691', '03742', '03898', '03907', '03915', '03975', 
# '04023', '04048', '04100', '04376', '04404', '04405', 
# '04409', '04458', '04470', '04479', '04537', '04641', 
# '04672', '04681', '04703', '04710', '05039', '05206', 
# '05325', '05359', '05447', '05535', '05717', '05780', 
# '05859', '05908', '06070', '06173', '06175', '06184', 
# '06207', '06233', '06256', '06261', '06264', '06286', 
# '06392', '06447', '06475', '06506', '06508', '06547', 
# '06723', '06758', '06841', '06861', '06924', '07204', 
# '07217', '07395', '07435', '07520', '07584', '07703', 
# '07775', '07815', '07837', '07975', '08022', '08049', 
# '08173', '08219', '08291', '08338', '08710', '08720', 
# '08721', '08738', '08742', '08975', '09061', '09075', 
# '09132', '09144', '09148', '09153', '09157', '09162', 
# '09176', '09177', '10001', '10004', '10006', '10008', 
# '10009', '10010', '10013', '10014', '10015', '10016', 
# '10017', '10018', '10019', '10020', '10021', '10022', 
# '10023', '10024', '10025', '10026', '10027', '10028', 
# '10029', '10030', '10031', '10032', '10033', '10034', 
# '10035', '10036', '10037', '10038', '10039', '10040', 
# '10041', '10042', '10043', '10044', '10045', '10046', 
# '10047', '10048', '10049', '10050', '10051', '10052', 
# '10053', '10054', '10055', '10056', '10058', '10059', 
# '10060', '10062', '10063', '10064')

cluster = c('10053')

myInputC = read.table("info.sort")
inputC = data.frame(cl = myInputC$V1, sigma = myInputC$V17, R200 = myInputC$V27, M200 = myInputC$V30)

for(cl in 1:length(cluster)) {
	dir = paste("clusters_nosocs/cluster.gals.sel.shiftgap.iter.", cluster[cl], sep="")
	myInput = read.table(dir)
	# myInput = read.table("sem_rotacao/fort.8")
	inputT = data.frame(ra = myInput$V1, dec = myInput$V2, zspec = myInput$V8, Vel = myInput$V10,flag = myInput$V15, r = myInput$V5)

	cluster_input = inputC[which(inputC$cl == as.integer(cluster[cl])),names(inputC) %in% c("sigma","R200","M200")]

	input <- inputT[which(inputT$flag == 0),names(inputT) %in% c("ra","dec","zspec","Vel","flag", "r")]
	if(nrow(input) <= 20)
		print(paste(cluster[cl], " Número objetos: ",nrow(input)), sep="")
	if(nrow(input) > 20) {
		fileName = paste('cluster_M', cluster[cl], '.txt', sep="")

		# file = paste("Aglomerado_rotacao M",ID, sep="")
		file = cluster[cl]

		dir.create(file)
		setwd(file)

		# Ordenar a partir da velocidade
		input = input[order(input$Vel),]

		write.table(input, file=fileName, sep=" ")

		########### ANÁLISE DE GAPS ###########

		jpeg('distribuição de velocidades.jpg')
		h = hist(input$Vel, main = "Distribuição de Velocidades", 
			     xlim = c(min(input$Vel),max(input$Vel)), breaks = 30)
		plot(h$mids, h$counts, type="s",
			 xlab=expression(paste(Velocidades,"  ","(","km","/",s^{-1},")")), ylab="N")

		xfit = seq(min(input$Vel), max(input$Vel), length=100)
		yfit = dnorm(xfit,mean=mean(input$Vel), sd = sd(input$Vel))
		yfit = yfit*diff(h$mids[1:2]*length(input$Vel))
		lines(xfit,yfit,col="blue",lwd=2)
		rug(input$Vel, ticksize=0.03, side=1, lwd=2, col="black")

		indexBCG = which.min(input$r) # posição do objeto mais brilhante

		# Plota no gráfico BCG
		abline(v = input$Vel[indexBCG],lty=5, lwd=2, col = "darkorchid4")
		legend((input$Vel[indexBCG]), max(h$counts),"BCG",
			    box.col = NA, bg = NA, adj = 0.2)

		legend("topright", inset=.05, legend=paste("N = ", nrow(input)), box.lty=0)

		# --- CALCULA O GAP ---
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
			print(paste("MEDIANA ", cluster[cl], " = ", medVel[1]))
		}

		# Insere no histograma os gaps > 2.25
		for(n in 1:length(gap_positions)) {
			aux = gap_positions[n]
			segments(x0=(input$Vel[aux]+((input$Vel[aux+1]-input$Vel[aux])/2)), 
				 y0=0.1, x1 = (input$Vel[aux]+((input$Vel[aux+1]-input$Vel[aux])/2)), 
		         y1 = 1.5, col="red", lwd=2.5)
		}

		title(main = file, sub = "",
		  cex.main = 2,   font.main= 1, col.main= "darkorchid4")

		dev.off() # fecha o histograma

		########### CALCULA A DIREÇÃO DO EIXO PRINCIPAL ###########

		# Encontra a posição do maior gap
		pos_gap = which.max(weiGap)

		input['Flag'] = 0

		# Valores à esquerda do maior gap recebe em Flag o valor -1 e à direita 1
		for(i in 1:N) {
			if(pos_gap >= i) 
				input$Flag[i] = -1
			else
				input$Flag[i] = 1
		} 

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

		# Insere os valores de x e y no dataframe
		input['Vx'] = x
		input['Vy'] = y

		jpeg('Eixo Principal.jpg') # cria arquivo de imagem

		eli = ellipse(cor(x,y),scale=c(sd(x),sd(y)),centre=c(mean(x),mean(y)),level=0.95, npoints=round(length(x)*0.95))
		aplot(eli[,1],eli[,2],type="l",asp=1,col='darkorchid1')

		# Esboça os valores de x e y (esquerda e direita do gap)
		points(x,y, pch=ifelse(input$Flag == -1, "-", "+"),
			  col=ifelse(input$Flag == -1, 'red', 'blue'), cex = 1.5)

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
		lines(dx, ypred, col="gold1", lwd=2)

		pnts = cbind(x,y)
		distAxis = apply(pnts, 1, function(x) dist_point_line(as.numeric(x[1:2]), slope = alpha, intercept = beta))
		draw.circle(center_elip[1],center_elip[2],radius=0.5, border = "darkorchid4")

		title(main = file, sub = "",
		  cex.main = 2,   font.main= 1, col.main= "darkorchid4")

		dev.off() # fecha elipse

		# Insere as distâncias na última coluna do dataframe
		input['DistA'] = distAxis

		########### TESTES ESTATÍSTICOS ###########
		rotacao = FALSE # considerando que não haja uma possível rotação

		con = file("Teste Espacial de duas amostras.txt")
		sink(con)

		# Todos os pontos do gráfico
		a1 = subset(input[c("Vx","Vy")], input$Flag == -1)
		a2 = subset(input[c("Vx","Vy")], input$Flag == 1)
		a1 = matrix(c(a1$Vx, a1$Vy), ncol = 2)
		a2 = matrix(c(a2$Vx, a2$Vy), ncol = 2)

		if(length(a1) > 0 && length(a2) > 0) {
		 	cat("------------------------ TODOS OS PONTOS ------------------------\n\n")
			cat("-------- TESTE DE CRAMER --------\n")
			cramer = cramer.test(a1,a2)
			print(cramer$p.value)
			
			cat("-------- TESTE DE HOTELLING --------\n")
			hotelling = hotelling.test(a1,a2)
			print(hotelling$pval)

			if((!is.na(cramer$p.value) && cramer$p.value < 0.05) || (!is.na(hotelling$pval) && hotelling$pval < 0.05))
				rotacao = TRUE
		}
		# Todos os pontos acima do eixo principal
		a1 = subset(input[c("Vx","Vy")], input$Flag == -1 & input$DistA >= 0)
		a2 = subset(input[c("Vx","Vy")], input$Flag == 1 & input$DistA >= 0)
		a1 = matrix(c(a1$Vx,a1$Vy), ncol=2)
		a2 = matrix(c(a2$Vx,a2$Vy), ncol=2)

		if(length(a1) > 0 && length(a2) > 0) {
			cat("\n\n------------------------ PONTOS ACIMA DO EIXO ------------------------\n\n")
			cat("-------- TESTE DE CRAMER --------\n")
			cramer = cramer.test(a1,a2)
			print(cramer$p.value)
			
			cat("-------- TESTE DE HOTELLING --------\n")
			hotelling = hotelling.test(a1,a2)
			print(hotelling$pval)
			
			if((!is.na(cramer$p.value) && cramer$p.value < 0.05) || (!is.na(hotelling$pval) && hotelling$pval < 0.05))
				rotacao = TRUE
		}

		# Todos os pontos abaixo do eixo principal
		a1 = subset(input[c("Vx","Vy")], input$Flag == -1 & input$DistA < 0)
		a2 = subset(input[c("Vx","Vy")], input$Flag == 1 & input$DistA < 0)
		a1 = matrix(c(a1$Vx,a1$Vy), ncol=2)
		a2 = matrix(c(a2$Vx,a2$Vy), ncol=2)

		if(length(a1) > 0 && length(a2) > 0) {
			cat("\n\n------------------------ PONTOS ABAIXO DO EIXO ------------------------\n\n")
			cat("-------- TESTE DE CRAMER --------\n")
			cramer = cramer.test(a1,a2)
			print(cramer$p.value)
			
			cat("-------- TESTE DE HOTELLING --------\n")
			hotelling = hotelling.test(a1,a2)
			print(hotelling$pval)

			if((!is.na(cramer$p.value) && cramer$p.value < 0.05) || (!is.na(hotelling$pval) && hotelling$pval < 0.05))
				rotacao = TRUE
		}

		sink() # fecha arquivo teste espacial de duas amostras

		########### CALCULAR O VALOR DE W ROTACIONAL ###########
		if(rotacao == TRUE) {
			print(paste(file,"com rotação"))
			# Cria um outro dataframe
			velDist = data.frame(Vel = input$Vel, DistA = abs(input$DistA), 
			 Flag = input$Flag, Vx = input$Vx, Vy = input$Vy)

			velDist[,1] = velDist[,1] - cluster_input$sigma
			velDist[,1] = velDist[,1] * 3.2408e-20 		# converte de km/s para mpc/s
			velDist = velDist[order(velDist[2]),]       # ordena a partir da distância

			# Soma as velocidades de 20 galáxias (+ e -)
			velNegat = 0
			velPosit = 0
			vel20 = velDist$Flag[1:20]	

			velNegat = c()
			velPosit = c()

			velNegat = velDist$Vel[which(vel20 == -1)]
			velPosit = velDist$Vel[which(vel20 == 1)]

			wiNegat = length(velNegat) + 1
			wiPosit = length(velPosit) + 1

			deltaV = c()
			deltaV[1] = median(velPosit) - median(velNegat)

			raio = c()
			raio[1] = velDist$DistA[20]

			w = c()
			w[1] = deltaV[1]/raio[1]
			j = 2

			for(i in 21:nrow(velDist)) {
				velocity = velDist$Vel[i]
		  		if(velDist$Flag[i] == -1) {
		  			velNegat[wiNegat] = velocity
		  			wiNegat = wiNegat + 1
		  		}
		  		else{
		  			velPosit[wiPosit] = velocity 
		  			wiPosit = wiPosit + 1
		  		}

		  		deltaV[j] = median(velPosit) - median(velNegat)
		  		raio[j] = velDist$Dist[i]
		  		w[j] = deltaV[j]/raio[j]
		  		j = j + 1
			}

				distance = numeric(0)
			for (i in 1:nrow(velDist)){
			    distance[i] = dist_2_points(center_elip, c(velDist[i,5],velDist[i,6]))
			}

			velDist['dist'] = distance

			subR = subset(velDist, velDist$dist <= 0.5)
				posBCGV = which.min(subR$r)
			
			write.table(subR, file="subr.txt", sep=" ")
			write.table(deltaV, file="deltaV.txt", sep=" ")

			velocityMax = omegaV(a_point, b_point, center_elip ,cluster_input$R200, cluster_input$M200 * 1e+14)

			jpeg('PerfilDeVelocidadeAngular.jpg')   # cria imagem para gráfico
			plot(raio,abs(w), type="l", xlab="Raio (Mpc)", ylab=expression(paste(omega~"(R)", " (rad/s)")),col="darkblue")
			abline(v = cluster_input$R200,lty=5, lwd=2, col="green")

			abline(h = velocityMax, lty=3, lwd=2, col="red")

			AF1 = approxfun(raio,abs(w))
			points(cluster_input$R200, AF1(cluster_input$R200), pch=20, col="black")  
			legend(cluster_input$R200, AF1(cluster_input$R200), signif(AF1(cluster_input$R200), digits = 6), box.col = NA, bg = NA, adj = 0.2)

			title(main = file, sub = "",
				cex.main = 2,   font.main= 1, col.main= "darkorchid4")
			dev.off() # fecha o arquivo de perfil de velocidade angular

			inputRotate = data.frame(distancia = velDist$dist, Vel = input$Vel)

			inputRotate[,2] = inputRotate[,2] * 3.2408e-20 	# converte de km/s para mpc/s
			inputRotate[,1] = inputRotate[,1] * 0.001	    # converte de kpc para mpc

			inputRotate = inputRotate[order(inputRotate[1]),]   # ordena a partir da distância
			
			jpeg('CurvaDeRotacao.jpg')   # cria imagem para gráfico
			plot(inputRotate$distancia,inputRotate$Vel, type="l", xlab="Distância (Mpc)", ylab="Velocidade (Mpc/s)",col="darkblue")

			title(main = file, sub = "",
				cex.main = 2,   font.main= 1, col.main= "darkorchid4")
			dev.off() # fecha o arquivo de perfil de velocidade angular
		}
		if(rotacao == FALSE)
			print(paste(file,"sem rotação"))

		setwd("../")
	}
}
