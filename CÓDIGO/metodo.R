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

myInput = read.table("selec20.dat")
# myInput = read.table("sem_rotacao/fort.8")
inputT = data.frame(ID = myInput$clusterID, ra = myInput$ra, dec = myInput$dec, zspec = myInput$zspec, flag = myInput$Flag, Vel = myInput$Vel)

# myInputRot = read.table(fileNameRot)

for(ID in 1:max(inputT$ID)){ 
	input <- inputT[which(inputT$ID == ID & inputT$flag == 0),names(inputT) %in% c("ID","ra","dec","zspec","Vel")]

	fileName = paste('cluster ', ID, '.txt', sep="")

	# file = paste("Aglomerado_rotacao M",ID, sep="")
	file = paste("Aglomerado ",ID, sep="")

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
	legend("top", inset=.05, legend=paste("N = ", nrow(input)), box.lty=0)

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
		print(paste("MEDIANA ", ID))
		medVel = input$Vel[which(input$Vel >= median(input$Vel))]
		gap_positions = medVel[1]
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
 		  col=ifelse(input$Flag == -1, 'firebrick1', 'chartreuse2'))

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

		velDist[,1] = velDist[,1] * 3.2408e-20 	# converte de km/s para mpc/s
		velDist = velDist[order(velDist[2]),]       # ordena a partir da distância

		# Soma as velocidades de 20 galáxias (+ e -)
		velNegat = 0
		velPosit = 0

		for(i in 1:20) {
			if(velDist$Flag[i] == -1) 
				velNegat = velDist$Vel[i] + velNegat
			else
				velPosit = velDist$Vel[i] + velPosit
		}
		deltaV = c()
		deltaV[1] = velPosit - velNegat

		raio = c()
		raio[1] = velDist$DistA[20]

		w = c()
		w[1] = deltaV[1]/raio[1]
		j = 2

		for(i in 21:nrow(velDist)) {
	  		if(velDist$Flag[i] == -1) {
	  			velNegat = velDist$Vel[i] + velNegat
	  		}
	  		else{
	  			velPosit = velDist$Vel[i] + velPosit
	  		}

	  		deltaV[j] = velPosit - velNegat
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
 		
 		jpeg('PerfilDeVelocidadeAngular.jpg')   # cria imagem para gráfico
  		plot(raio,abs(w), type="l", xlab="Raio (Mpc)", ylab=expression(paste(omega~"(R)", " (Mpc/s)")),col="darkblue")

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
	else {
		print(paste(file,"sem rotação"))
	}

	setwd("../")
}