rotate <- function(text) {
	library(stats)
	library(astro)
	library(ellipse)
	library(cosmoFns)
	library(cramer)
	library(Hotelling)
	library(perm)

# --- LEITURA DE DADOS ---
	myInput = read.table(text)
	maxID = 20
	#maxID = max(myInput$clusterID) 
	for(ID in 1:maxID) {
		# Insere no dataframe os objetos contidos em cada aglomerado (flag)
		input = data.frame(clusterID = myInput$clusterID, ra = myInput$ra, dec = myInput$dec,
		                   r = myInput$r, zspec = myInput$zspec, Vel = myInput$Vel,
		                   Flag = myInput$Flag, lgm_tot_p50 = myInput$lgm_tot_p50)
		input = subset(input, (myInput$clusterID == ID & myInput$Flag == 0))
		#input = subset(input, lgm_tot_p50 < 10.25) #massa

		# Ordenar a partir da velocidade
		input = input[order(input$Vel),]

# --- CONSTRUÇÃO DO HISTOGRAMA DE DISTRIBUIÇÃO DE VELOCIDADES ---
		
		# Cria pasta do aglomerado ID
		file = paste("Aglomerado MASS_I",ID)
		dir.create(file)
		setwd(file)

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
		abline(v = input$Vel[indexBCG],lty=5, lwd=2)
		legend((input$Vel[indexBCG]-500), median(h$counts),"BCG",
			    box.col = "white", bg = "white", adj = 0.2)

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

		# Insere no histograma os gaps > 2.25
		for(n in 1:length(gap_positions)) {
			aux = gap_positions[n]
			segments(x0=(input$Vel[aux]+((input$Vel[aux+1]-input$Vel[aux])/2)), 
				 y0=0.1, x1 = (input$Vel[aux]+((input$Vel[aux+1]-input$Vel[aux])/2)), 
                 y1 = 1.5, col="red", lwd=2.5)
		}
		dev.off() # fecha o histograma

# --- CALCULA A DIREÇÃO DO EIXO PRINCIPAL ---

		# Encontra a posição do maior gap
		pos_gap = which.max(weiGap)

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

# CALCULA A DIREÇÃO DO EIXO PRINCIPAL ---
		jpeg('Eixo Principal.jpg') # cria arquivo de imagem

		eli = ellipse(cor(x,y),scale=c(sd(x),sd(y)),centre=c(mean(x),mean(y)),level=0.95, npoints=round(length(x)*0.95))
		aplot(eli[,1],eli[,2],type="l",asp=1,col='darkorchid1')
 	
	 	# Esboça os valores de x e y (esquerda e direita do gap)
	 	points(x,y, pch=ifelse(input$Flag == -1, "-", "+"),
	 		  col=ifelse(input$Flag == -1, 'firebrick1', 'chartreuse2'))

	 	points(input$Vx[indexBCG],input$Vy[indexBCG],pch = 5,col="blue",cex=1.5 )
	 	legend(input$Vx[indexBCG], input$Vy[indexBCG], "BCG",box.col = "transparent",bg = "transparent", adj = 0.2)

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

 		# Insere as distâncias na última coluna do dataframe
	 	input['DistA'] = distAxis
	 	write.table(distAxis, file="dist.txt", sep=" ")

# --- APLICAR TESTE DE SIMETRIA ---
	con = file("Teste de Simetria.txt")
	sink(con)
	cat("-------- TESTE DE SIMETRIA À DIREITA DO GAP --------\n\n")
	cat("------------------------ TODOS OS PONTOS ------------------------\n\n")

	# Todos os pontos do gráfico à direita do gap
	allDpermD = subset(input, input$Flag == 1 & input$DistA >  0)
	allDpermE = subset(input, input$Flag == 1 & input$DistA <  0)
	permt = permTS(abs(allDpermD$DistA),abs(allDpermE$DistA))
	print(permt$p.value)

	# Todos os pontos do gráfico à esquerda do gap
	allEpermD = subset(input, input$Flag == -1 & input$DistA >  0)
	allEpermE = subset(input, input$Flag == -1 & input$DistA <  0)
	permt = permTS(abs(allEpermD$DistA),abs(allEpermE$DistA))
	print(permt$p.value)

	cat("\n\n------------------------ PONTOS ACIMA DO EIXO ------------------------\n\n")
	# Todos os pontos acima do eixo principal
	abovepermD = subset(input, input$Flag == 1 & input$DistA >  0)
	abovepermE = subset(input, input$Flag == -1 & input$DistA >  0)
	permt = permTS(abs(abovepermD$DistA),abs(abovepermE$DistA))
	print(permt$p.value)

	cat("\n\n------------------------ PONTOS ABAIXO DO EIXO ------------------------\n\n")
	# Todos os pontos abaixo do eixo principal
	belowpermD = subset(input, input$Flag == 1 & input$DistA <  0)
	belowpermE = subset(input, input$Flag == -1 & input$DistA <  0)
	permt = permTS(abs(belowpermD$DistA),abs(belowpermE$DistA))
	print(permt$p.value)

	sink() # fecha arquivo teste de simetria

# --- APLICAR OS TESTES DE CRAMER E HOTELLING ---
	 	con = file("Teste Espacial de duas amostras.txt")
	 	sink(con)

	 	# Todos os pontos do gráfico
	 	a1 = subset(input[c("Vx","Vy")], input$Flag == -1)
	 	a2 = subset(input[c("Vx","Vy")], input$Flag == 1)
	 	a1 = matrix(c(a1$Vx, a1$Vy), ncol = 2)
	 	a2 = matrix(c(a2$Vx, a2$Vy), ncol = 2)
	 	cat("------------------------ TODOS OS PONTOS ------------------------\n\n")
		cat("-------- TESTE DE CRAMER --------\n")
		cramer = cramer.test(a1,a2)
		print(cramer$p.value)
		
		cat("-------- TESTE DE HOTELLING --------\n")
		hotelling = hotelling.test(a1,a2)
		print(hotelling$pval)

		# Todos os pontos acima do eixo principal
		a1 = subset(input[c("Vx","Vy")], input$Flag == -1 & input$DistA >= 0)
		a2 = subset(input[c("Vx","Vy")], input$Flag == 1 & input$DistA >= 0)
		a1 = matrix(c(a1$Vx,a1$Vy), ncol=2)
		a2 = matrix(c(a2$Vx,a2$Vy), ncol=2)
		
		cat("\n\n------------------------ PONTOS ACIMA DO EIXO ------------------------\n\n")
		cat("-------- TESTE DE CRAMER --------\n")
		cramer = cramer.test(a1,a2)
		print(cramer$p.value)
		
		cat("-------- TESTE DE HOTELLING --------\n")
		hotelling = hotelling.test(a1,a2)
		print(hotelling$pval)
		
		# Todos os pontos abaixo do eixo principal
		a1 = subset(input[c("Vx","Vy")], input$Flag == -1 & input$DistA < 0)
		a2 = subset(input[c("Vx","Vy")], input$Flag == 1 & input$DistA < 0)
		a1 = matrix(c(a1$Vx,a1$Vy), ncol=2)
		a2 = matrix(c(a2$Vx,a2$Vy), ncol=2)
		
		cat("\n\n------------------------ PONTOS ABAIXO DO EIXO ------------------------\n\n")
		cat("-------- TESTE DE CRAMER --------\n")
		cramer = cramer.test(a1,a2)
		print(cramer$p.value)
		
		cat("-------- TESTE DE HOTELLING --------\n")
		hotelling = hotelling.test(a1,a2)
		print(hotelling$pval)

		sink() # fecha arquivo teste espacial de duas amostras

# --- CALCULAR O VALOR DE W ROTACIONAL ---

		# Cria um outro dataframe
		velDist = data.frame(Vel = input$Vel, DistA = abs(input$DistA), 
			 Flag = input$Flag, r = input$r, Vx = input$Vx, Vy = input$Vy)

		velDist[,1] = velDist[,1] * 3.24078e-20 	# converte de km/s para mpc
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
  		info = read.table("../info.txt")
  		R200 = info$V18[ID]
  		distance = numeric(0)
		for (i in 1:nrow(velDist)){
		    distance[i] = dist_2_points(center_elip, c(velDist[i,5],velDist[i,6]))
		}
		velDist['dist'] = distance

		subR = subset(velDist, velDist$dist <= 0.5)
  		posBCGV = which.min(subR$r)
  		
  		points(subR$Vx[posBCGV],subR$Vy[posBCGV],pch = 5,col="red",cex=1.5 )
 		legend(subR$Vx[posBCGV],subR$Vy[posBCGV], "BCGV",box.col = "transparent",bg = "transparent", adj = 0.2)
 		dev.off() # fecha arquivo imagem eixo principal
		
		write.table(subR, file="subr.txt", sep=" ")
 		jpeg('PerfilDeVelocidadeAngular.jpg')   # cria imagem para gráfico
  		plot(raio,abs(w), type="l", xlab="Raio(Mpc)", ylab=expression(omega~"(R)"),col="darkblue")
  		abline(v = R200,lty=5, lwd=2, col="green")
  		
  		dev.off() # fecha o arquivo de perfil de velocidade angular
 		

		setwd("../")
	} # fecha laço ID do aglomerado

	
}


# função distâncias entre dois pontos
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
