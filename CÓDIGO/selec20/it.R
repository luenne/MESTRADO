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


myInput = read.table('selec20.dat')
inputT = data.frame(clusterID = myInput$clusterID, ra = myInput$ra, dec = myInput$dec,
		                   r = myInput$r, zspec = myInput$zspec, Vel = myInput$Vel,
		                   Flag = myInput$Flag)

sample = c()
rot = matrix(nrow = 1000, ncol = max(inputT$clusterID))
for(it in 1:1000) {
	for(ID in 1:max(inputT$clusterID)) {

		input <- subset(inputT, (inputT$clusterID == ID & inputT$Flag == 0))

		# Ordenar a partir da velocidade
		input = input[order(input$Vel),]

		########### ANÁLISE DE GAPS ###########

		indexBCG = which.min(input$r) # posição do objeto mais brilhante

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
			#print(paste("MEDIANA ", ID, " = ", medVel[1]))
		}

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

		########### TESTES ESTATÍSTICOS ###########
		rotacao = FALSE # considerando que não haja uma possível rotação

		# Todos os pontos do gráfico
		a1 = subset(input[c("Vx","Vy")], input$Flag == -1)
		a2 = subset(input[c("Vx","Vy")], input$Flag == 1)
		a1 = matrix(c(a1$Vx, a1$Vy), ncol = 2)
		a2 = matrix(c(a2$Vx, a2$Vy), ncol = 2)
		
		cramer = cramer.test(a1,a2)
		hotelling = hotelling.test(a1,a2)

		if(length(a1) > 0 && length(a2) > 0) {

			if((!is.na(cramer$p.value) && cramer$p.value < 0.05) || (!is.na(hotelling$pval) && hotelling$pval < 0.05)){
				rotacao = TRUE
			}
		}
		# Todos os pontos acima do eixo principal
		a1 = subset(input[c("Vx","Vy")], input$Flag == -1 & input$DistA >= 0)
		a2 = subset(input[c("Vx","Vy")], input$Flag == 1 & input$DistA >= 0)
		a1 = matrix(c(a1$Vx,a1$Vy), ncol=2)
		a2 = matrix(c(a2$Vx,a2$Vy), ncol=2)

		cramer = cramer.test(a1,a2)
		hotelling = hotelling.test(a1,a2)

		if(length(a1) > 0 && length(a2) > 0) {
			
			if((!is.na(cramer$p.value) && cramer$p.value < 0.05) || (!is.na(hotelling$pval) && hotelling$pval < 0.05)){
				rotacao = TRUE
			}
		}

		# Todos os pontos abaixo do eixo principal
		a1 = subset(input[c("Vx","Vy")], input$Flag == -1 & input$DistA < 0)
		a2 = subset(input[c("Vx","Vy")], input$Flag == 1 & input$DistA < 0)
		a1 = matrix(c(a1$Vx,a1$Vy), ncol=2)
		a2 = matrix(c(a2$Vx,a2$Vy), ncol=2)

		cramer = cramer.test(a1,a2)
		hotelling = hotelling.test(a1,a2)

		if(length(a1) > 0 && length(a2) > 0) {

			if((!is.na(cramer$p.value) && cramer$p.value < 0.05) || (!is.na(hotelling$pval) && hotelling$pval < 0.05)){
				rotacao = TRUE
			}
		}

		sample[ID] = rotacao
	}
	print(it)
	rot[it, ] = sample 

}
