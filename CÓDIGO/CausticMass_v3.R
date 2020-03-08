
biweightScale=function(x, tuningConstant=9){
  Median=median(x)
  MAD=mad(x, constant = 1)
  len=length(x)
  diffModuli=vector(length = len)
  for (i in 1:len) {
    diffModuli[i]=abs(x[i]-Median)
  }
  uValues=vector(length = len)
  for (i in 1:len) {
    uValues[i]=abs(x[i]-Median)/(tuningConstant*MAD)
  }
  top=bottom=valCount=0
  for (i in 1:len) {
    if(abs(uValues[i]) <= 1){
      u2Term=1-uValues[i]^2
      u4Term=u2Term^4
      top=top+diffModuli[i]^2*u4Term
      bottom=bottom+(u2Term*(1-5*uValues[i]^2))
      valCount=valCount+1
    }
  }
  top=sqrt(top)
  bottom=abs(bottom)
  SBI=sqrt(valCount)*(top/bottom)
  return(SBI)
}

shiftgapper=function(dproj, vlos, npbin=25, gapper=5, plot=T){
  require(OneR)
  require(robustbase)
  
  gap_prev = 2000
  nsize=ngal=length(dproj)
  names(vlos)=1:nsize
  vlos.=vlos[order(dproj)]; dproj.=sort(dproj)
  nbins=round(ceiling(nsize/npbin))
  nbins=ifelse(nbins > 1, nbins, 2)
  bindata=OneR::bin(dproj., nbins, 1:nbins, 'content')
  datanew=datafinal=c()
  
  for (i in 1:nbins) {
    #print(i)
    databin=vlos.[bindata == i]
    nsize=length(databin)
    datasize = nsize-1
    if(nsize > 5){
      while(nsize - datasize > 0 & datasize >= 5){
        nsize=length(databin)
        databinsort=sort(databin)
        f=switch(gapper,
                 '1'=(databinsort[nsize-ceiling(nsize/4)]-databinsort[ceiling(nsize/4)])/1.349,
                 '2'=(quantile(databinsort, .75)-quantile(databinsort, .25))/1.349,
                 '3'=robustbase::s_Qn(databinsort),
                 '4'=biweightScale(databinsort),
                 '5'=1000)
        #gap=f/1.349
        gap=f
        if (gap < 500) break;
        databelow=databinsort[databinsort <= 0]
        nbelow=length(databelow)
        if(nbelow != 0){
          gapbelow=databelow[2:nbelow]-databelow[1:(nbelow-1)]
        }else{
          gapbelow=-1; databelow=logical(0)
        }
        if(max(gapbelow, na.rm = T) >= gap){
          vgapbelow=which(gapbelow >= gap)[1]
        }else{
          vgapbelow=0
        }
        datanew=c(datanew, databelow[(vgapbelow+1):nbelow])
        dataabove=databinsort[databinsort > 0]
        nabove=length(dataabove)
        if(nabove != 0){
          gapabove=dataabove[2:nabove]-dataabove[1:(nabove-1)]
        }else{
          gapabove=-1; dataabove=logical(0)
        }
        if(max(gapabove, na.rm = T) >= gap){
          vgapabove=which(gapabove >= gap)[1]
        }else{
          vgapabove = nabove
        }
        datanew=c(datanew, dataabove[1:vgapabove])
        databin = datanew
        datasize=length(datanew)
        datanew=c()
      }
      gap_prev=ifelse(gap >= 500, gap, 500)
    }
    datafinal=c(datafinal, databin)
  }
  namesdata=as.numeric(names(datafinal))
  gclass=rep(1, ngal); gclass[namesdata]=0
  if(plot){
    graphics::plot(dproj, vlos, pch=20, col=gclass+1)
  }
  return(gclass)
}

gaussian_kernel=function(dproj, vlos, r200, norm=100, scale=10, xmax=6, ymax=5e3, by=0.05, plot=F){
  # Uses a 2D gaussian kernel to estimate the density of the phase space.
  # As of now, the maximum radius extends to 6Mpc and the maximum velocity allowed is 5000km/s
  # The "q" parameter is termed "scale" here which we have set to 10 as default, but can go as high as 50.
  # "normalization" is simply H0
  # "x/yres" can be any value, but are recommended to be above 150
  # "adj" is a custom value and changes the size of uniform filters when used (not normally needed)
  require(gplots)
  require(magicaxis)
  require(spatstat)
  
  if(any(dproj >= xmax)){
    stop(paste('Error: Please either increase your xmax value or trim your sample to be x < ', xmax))
  } 
  if(any(abs(vlos) >= ymax)){
    stop(paste('Error: Please either increase your ymax value or trim your sample to be y < ', ymax))
  } 
  normscale=norm*scale
  #xvalues=rep(dproj, 2)
  #yvalues=c(vlos, -vlos)
  xvalues=dproj; yvalues=vlos
  yvalues=yvalues/normscale
  x_range=seq(0, xmax-by, by=by)
  xres=length(x_range)
  y_range=seq(-ymax/normscale, ymax/normscale-by, by=by)*normscale
  yres=length(y_range)
  x_scale = (xvalues/xmax)*xres
  y_scale=((yvalues*normscale+ymax)/(2*ymax))*length(y_range)
  if(length(which(xvalues < r200)) < 7){ # NEW
    r200=unique(sort(xvalues))[7]
    #message('Too few elements inside r200, extending radius limit')
    warning('Too few elements inside r200, extending radius limit', immediate. = T, call. = F)
  }
  bwx=biweightScale(x_scale[xvalues<r200], 9)
  bwy=biweightScale(y_scale[xvalues<r200], 9)
  ksize=(4/(3*length(xvalues)))^(1/5)*sqrt((bwx^2+bwy^2)/2)
  ymx=ymax/normscale
  h=gplots::hist2d(c(xvalues,0,xmax,xmax,0), c(yvalues,-ymx,-ymx,ymx,ymx), c(xres,yres), show = F)
  h$counts[1,1]=h$counts[xres,1]=h$counts[xres,yres]=h$counts[1,yres]=0
  img=spatstat::blur(spatstat::as.im(h$counts), sigma = ksize, normalise = T)$v
  #img=MASS::kde2d(xvalues, yvalues*normscale, n=c(xres, yres), lims = c(range(x_range),range(y_range)))$z
  #img=fields::image.smooth(h$counts)$z
  #img=as.matrix(imager::isoblur(imager::as.cimg(h$counts), sigma = ksize))
  #img=as.matrix(imager::deriche(imager::as.cimg(h$counts), sigma=ksize_x, order=2, axis="y"))
  #img=applyFilter(h$counts, kernel = convKernel(sigma = ksize_x, k = 'gaussian'))
  #img=raster::as.matrix(spatialEco::raster.gaussian.smooth(raster::raster(c$imgr), 4))
  if(plot){
    magicaxis::magimage(x_range, y_range, img, asp = NA, xlab=expression(R[proj]~(Mpc)), 
                        ylab=expression(v[proj]~(km~s^{-1})))
  }
  lout=list(x=x_range, y=y_range, z=img, sigma=ksize)
  return(lout)
}

restrict_gradient=function(pastA,newA,pastr,newr){
  # It is necessary to restrict the gradient the caustic can change at in order to be physical
  gradu = 0.5; gradd = 2.0
  if(pastA <= newA){
    if((log(newA)-log(pastA))/(log(newr)-log(pastr)) > gradu & pastA != 0){
      dr = log(newr)-log(pastr)
      return(exp(log(pastA) + gradu*dr))
    }else{
      return(newA)
    }
  }
  if(pastA > newA){
    if((log(newA)-log(pastA))/(log(newr)-log(pastr)) < -gradd & pastA != 0){
      dr = log(newr)-log(pastr)
      return(exp(log(pastA) - gradd*dr))
    }else{
      return(newA)
    }
  }
}

findcontours=function(Zi, ri, vi, r200, vvar, ri.max=4, nlevels=200, plot=T, verbose=T){
  #cl=contourLines(ri, vi, Zi, nlevels = nlevels)
  levels=10^(seq(0, log10(min(Zi[Zi > 0]/5)), length.out = nlevels))
  cl=contourLines(ri, vi, Zi, levels = levels)
  contours=list()
  for (i in 1:length(cl)) {
    #only consider contours that are "full" and don't loop back only in positive or negative space
    if(max(cl[[i]]$x) >= ri.max & min(cl[[i]]$x) <= 0 & max(cl[[i]]$y) > 0 & 
       min(cl[[i]]$y) < 0 & max(abs(cl[[i]]$y)) < 4500){
      #print(i)
      xcont_u=cl[[i]]$x[cl[[i]]$y > 0] #find positive/negative contours
      ycont_u=cl[[i]]$y[cl[[i]]$y > 0] 
      xcont_d=cl[[i]]$x[cl[[i]]$y < 0]
      ycont_d=cl[[i]]$y[cl[[i]]$y < 0]
      
      y_u=y_d=y_f=vector(length = length(ri)) #initialize positive, negative, and final arrays
      for (k in 1:length(ri)) { #loop over r grid values
        #match contour grid to r grid (nearest neighbor interpolate)
        wu=which(xcont_u > (ri[k]-0.05) & xcont_u < (ri[k]+0.05))
        y_u[k]=ifelse(length(wu) != 0, max(ycont_u[wu]), 0)
        wd=which(xcont_d > (ri[k]-0.05) & xcont_d < (ri[k]+0.05))
        y_d[k]=ifelse(length(wd) != 0, max(ycont_d[wd]), 0)
        
        if(k != 1){ #apply gradient restriction for positive and negative cases.
          y_u[k]=restrict_gradient(abs(y_u[k-1]), abs(y_u[k]), ri[k-1], ri[k])
          y_d[k]=restrict_gradient(abs(y_d[k-1]), abs(y_d[k]), ri[k-1], ri[k])
        }
        y_f[k] = min(y_u[k], abs(y_d[k])) #take minimum value of positive and negative arrays
      }
      contours=c(contours, list(y_f))
    }
  }
  if(length(contours) == 0) stop('The contours do not expand to the radial limit')
  
  #now I need to do the average calculation in Diaferio 99
  #because an integral is involved, I don't want to do this for all contours.
  #instead I select the 25% around the preliminary closest average and do
  #the full calculation for them
  avg_contours=unlist(lapply(contours, function(w) mean(w[ri < r200]^2)))
  avg_cont_diff=(avg_contours - 4*vvar)^2 #prelim diff calc
  i_sort=order(avg_cont_diff) #sort indices based on prelim diff
  i_sort_small=head(i_sort, as.integer(length(i_sort)/4))
  tot_avg=rep(0, length(i_sort_small))
  for (i in 1:length(i_sort_small)) {
    Ar = contours[[i]]
    lessr200 = which(ri <= r200)
    useri = ri[lessr200]
    Ar = Ar[lessr200]
    phir = rep(0, length(useri))
    for (j in 1:length(useri)) {
      philimit = abs(Ar[j])
      phir[j] = sum(Zi[j,][which((vi<philimit) & (vi>-philimit))])
    }
    tot_avg[i] = sum(Ar^2 * phir) / sum(phir)
  }
  final_contour = contours[[order((tot_avg - 4*vvar)^2)[1]]]
  
  if(plot){
    magimage(ri, vi, Zi, asp = NA, xlab=expression(R[proj]~(Mpc)), 
             ylab=expression(v[proj]~(km~s^{-1})), col.ticks = 'white')
    ll=lapply(contours, function(x) lines(ri, x, col='#FFBE74',lty=3))
    ll=lapply(contours, function(x) lines(ri, -x, col='#FFBE74',lty=3))
    lines(ri, final_contour, col='#F13005', lwd=2)
    lines(ri, -final_contour, col='#F13005', lwd=2)
  }
  if(verbose) message('complete')
  lout=list(caustic=final_contour, contours=contours)
  return(lout)
}

NFWfit=function(rii, Ar, halo_srad, ri_full){
  
  min_func=function(x, d0, halo_srad){
    sqrt(2*4*pi*4.5e-48*d0*halo_srad^2*log(1+x/halo_srad)/(x/halo_srad))*3.08e19
  }
  
  options(warn = -1)
  out=tryCatch(nls(y ~ min_func(x, d0, halo_srad=halo_srad), data=data.frame(x=rii, y=Ar), 
                   start = list(d0=1e14)), error=function(e) NA)
  options(warn = 0)
  if(anyNA(out)){
    out=nls(y ~ min_func(x, d0, halo_srad=halo_srad), data=data.frame(x=rii, y=Ar), start = list(d0=1e10))
  }
  halo_scale_density=coef(out)
  halo_scale_density_e=tryCatch(sqrt(halo_scale_density), error=function(e) 1e14)
  if(ri_full[1] == 0) ri_full[1]=ri_full[2]
  profile=min_func(ri_full, halo_scale_density, halo_srad)
  return(profile)
}

findvdisp=function(r, v, r200, maxv){
  # Use biweight sigma clipping Scale estimator for the velocity dispersion
  v_cut=v[which(r < r200 & abs(v) < maxv)]
  gal_vdisp=tryCatch(biweightScale(v_cut[is.finite(v_cut)], 9), error=function(e) NA)
  if(is.na(gal_vdisp)) gal_vdisp=sd(v_cut, na.rm = T)
  return(gal_vdisp)
}

findsurface=function(data, ri, vi, Zi, r200=2, maxv=5e3, memberflags=NA, halo_scale_radius=NA,
                      halo_scale_radius_e=0.01, halo_vdisp=NA, bin=NA, beta=NA, mirror=T, q=10, 
                      Hz=100, edge_perc=0.1, rimax=4, edge_int_remove=F, plot=T, verbose=T){
  require(stats)
  kappaguess=max(Zi)
  
  if(is.na(halo_scale_radius)){
    halo_scale_radius=r200/5
  }else{
    halo_scale_radius=halo_scale_radius
    halo_scale_radius_e=halo_scale_radius_e
  }
  
  # Calculate velocity dispersion with either members, fed value, or estimate using 3.5sigma clipping ----
  if(!is.na(memberflags)){
    vvarcal=data[,2][which(memberflags == 1)]
    gal_vdisp=tryCatch(biweightScale(vvarcal[vvarcal], 9), error=function(e) NA)
    if(!is.na(gal_vdisp) & verbose) print('O ya! membership calculation!')
    if(is.na(gal_vdisp)) gal_vdisp=sd(vvarcal)
    vvar=gal_vdisp^2
  }else{
    if(!is.na(halo_vdisp)){
      gal_vdisp=halo_vdisp
      vvar=gal_vdisp^2
    }else{
      gal_vdisp=tryCatch(findvdisp(data[,1], data[,2], r200, maxv), error=function(e) NA)
      if(is.na(gal_vdisp)) gal_vdisp=sd(data[,2][data[,1] < r200 & abs(data[,2]) < maxv])
      vvar=gal_vdisp^2
    }
  }
  
  # find contours ----------------------------------------------------------------------------------------
  fcont=findcontours(Zi, ri, vi, r200, vvar, ri.max = rimax, plot = plot, verbose = verbose)
  Ar_finalD=fcont$caustic
  contours=fcont$contours
  
  #print(Ar_finalD)
  
  # remove outliers from edge calculation ----------------------------------------------------------------
  # to implement!
  
  # Identify sharp phase-space edge ----------------------------------------------------------------------
  numbins = 6
  perc_top = edge_perc #what percent of top velocity galaxies per/bin used to identify surface
  numrval=length(which(data[,1] < r200)) #number of galaxies less than r200
  #if(numrval <= numbins) numbins=round(numrval/2)
  if(numrval == 1) numrval=length(which(data[,1] < 2*r200)) #NEW
  if(numrval <= numbins) numbins=numrval-1 #NEW
  size_bin=ceiling(numrval/numbins) #how many galaxies are in each bin
  rsort=sort(data[,1]) #sort r positions
  
  if(mirror){
    vsort=abs(data[,2][order(data[,1])]) #sort absolute value of velocities by r position
  }else{
    vsort=data[,2][order(data[,1])] #same as above but not abs
  }
  
  mid_rbin=avgmax=avgmin=mincomp=c()
  for (nn in 1:numbins) {
    #rng=((nn-1)*size_bin):((nn*size_bin)-1)
    rng=(size_bin*(nn-1)+1):(nn*size_bin)
    #print(rng)
    vbin=vsort[rng] #pick velocities in bin # nn
    if(length(vbin) == 0 & nn >= 4) break;
    rbin=rsort[rng] #pick radii in bin # nn
    #sort by velocity -> flip array from max-min -> take first edge_perc values where v>0
    vemax=head(sort(vbin, decreasing = T), ceiling(length(vbin[vbin>0])*perc_top))
    #sort by velocity -> take first edge_perc values where v<0
    vemin=head(sort(vbin), ceiling(length(vbin[vbin<0])*perc_top))
    avgmax=c(avgmax, mean(vemax, na.rm = T)) #add average of top edge_perc velocities to max array
    avgmin=c(avgmin, mean(vemin, na.rm = T)) #same as above but min array
    #take the minimum of either the above || below zero caustic
    if(is.nan(rev(avgmax)[1])) break;
    #if(is.nan(mean(vemax)) | is.na(mean(vemax))) break;
    #if no negative velocities (aka, mirrored), else take the minimum extreme
    #mincomp=ifelse(min(vbin) >= 0, c(mincomp, avgmax[nn]), c(mincomp, avgmin[nn]))
    if(min(vbin, na.rm = T) >= 0){
      mincomp=c(mincomp, avgmax[nn])
    }else{
      mincomp=c(mincomp, avgmin[nn])
    }
    mid_rbin=c(mid_rbin, median(rbin, na.rm = T)) #take median rvalue of bin
  }
  lcont=length(contours)
  chi=vector(length = lcont)
  #loop through contours and find squared difference with edge extreme
  for (nn in 1:lcont) {
    fint=stats::approxfun(ri[ri < r200], contours[[nn]][ri < r200]) #interpolate contour
    Ar_comp=fint(mid_rbin[mid_rbin < max(ri[ri < r200])]) #interpolated contour
    if(sum(Ar_comp) == 0) break;
    chi[nn]=median(abs(Ar_comp-mincomp[mid_rbin < max(ri[ri < r200])])) #measure squared distance
  }
  chi[chi == 0]=max(chi)
  #find level with min chi value
  Ar_finalE=tryCatch(contours[[which(chi == min(chi, na.rm = T))[1]]], error=function(e) rep(0, length(ri)))
  
  # fit an NFW to the resulting caustic profile ----------------------------------------------------------
  if(is.na(beta)){
    beta=rep(0.2, length(ri))
  }else{
    beta=beta
  }
  gb=(3-2*beta)/(1-beta)
  fitting_radii=which((ri>=r200/3) & (ri<=r200))
  rii=ri[fitting_radii]
  ArD=Ar_finalD[fitting_radii]*sqrt(gb[fitting_radii])
  ArE=Ar_finalE[fitting_radii]*sqrt(gb[fitting_radii])
  halo_srad=halo_scale_radius
  ri_full=ri
  vesc_fit=NFWfit(rii, ArD, halo_srad, ri_full)/sqrt(gb)
  vesc_fit_e=NFWfit(rii, ArE, halo_srad, ri_full)/sqrt(gb)
  
  # Output galaxy membership -----------------------------------------------------------------------------
  kpc2km = 3.09e16
  fitfunc=function(x, a, b) sqrt(2*4*pi*6.67e-20*a*(b*kpc2km)^2*log(1+x/b)/(x/b))
  reg=tryCatch(nls(y ~ fitfunc(x, a, b), data=data.frame(x=ri, y=Ar_finalD), start = list(a=5e14, b=1)),
               error=function(e) F)
  if(!reg){
    reg=tryCatch(nls(y ~ fitfunc(x, a, b=30), data=data.frame(x=ri, y=Ar_finalD), start = list(a=5e14)),
                 error=function(e) F)
  }
  
  memflag=vector(length = nrow(data))
  fcomp=stats::approxfun(ri, vesc_fit)
  for (k in 1:nrow(data)) {
    vcompare = fcomp(data[k,1])
    if(abs(vcompare) > abs(data[k,2])) memflag[k]=1
  }
  
  # plotting ---------------------------------------------------------------------------------------------
  if(plot){
    points(data[,1], data[,2], pch=4, col=3, cex=0.8)
    lines(ri, Ar_finalE, col='#4DC4CA', lwd=1, lty=2)
    lines(ri, -Ar_finalE, col='#4DC4CA', lwd=1, lty=2)
    lines(ri_full, vesc_fit_e, col='#E7D4C5', lwd=2, lty=2)
    lines(ri_full, -vesc_fit_e, col='#E7D4C5', lwd=2, lty=2)
    lines(ri_full, vesc_fit, col='#D0F1DC', lwd=2, lty=2)
    lines(ri_full, -vesc_fit, col='#D0F1DC', lwd=2, lty=2)
    points(data[,1][memflag == 0], data[,2][memflag == 0], pch=4, col='#EEEBE4', cex=0.8)
  }
  
  lout=list(caustic_profile=Ar_finalD, caustic_fit=vesc_fit, caustic_edge=abs(Ar_finalE), 
            caustic_fit_edge=vesc_fit_e, gal_vdisp=gal_vdisp, memflag=memflag)
  return(lout)
}

MassCalc=function(ri, A, vdisp, clus_z, r200=NA, conc1=NA, beta=0.25, fbr=NA, H0=100, Om=0.25, Ol=0.75){
  require(pracma)
  
  G = 6.67E-11
  solmass = 1.98892e30
  crit = 2.7745946e11*(H0/100.0)^2*(Om*(1+clus_z)^3 + Ol)
  r2 = ri[ri>=0]
  A2 = A[ri>=0]
  kmMpc = 3.08568025e19
  
  if(is.na(conc1)){
    conc = 5 + rnorm(1, 0, 2)
    if(conc < 0){
      conc = 5
    }
  }else{
    conc=conc1
  }
  beta=0.5*(ri/(ri+r200/conc))
  g_b=(3-2.0*beta)/(1-beta)
  
  sumtot=rep(0, length(A2))
  if(is.na(fbr)){
    f_beta = 0.5*((r2/r200*conc)^2)/((1+(r2/r200*conc))^2*log(1+((r2/r200*conc))))*g_b
    f_beta[1] = 0
    for (i in 1:(length(A2)-1)) {
      sumtot[i+1]=pracma::trapz(r2[2:(i+1)]*kmMpc*1000,f_beta[2:(i+1)]*(A2[2:(i+1)]*1000)^2)
      #sumtot[i+1]=Bolstad2::sintegral(r2[1:(i+1)]*kmMpc*1000,f_beta[1:(i+1)]*(A2[1:(i+1)]*1000)^2)$int
      #sumtot[i+1]=MESS::auc(r2[1:(i+1)]*kmMpc*1000,f_beta[1:(i+1)]*(A2[1:(i+1)]*1000)^2, type = 'spline')
    }
  }else{
    if(length(fbr) == 1){
      f_beta=rep(fbr, length(A2))
    }else{
      f_beta=fbr
    }
    f_beta[1] = 0
    for (i in 1:(length(A2)-1)) {
      sumtot[i+1]=pracma::trapz(r2[2:(i+1)]*kmMpc*1000, f_beta[2:(i+1)]*(A2[2:(i+1)]*1000)^2)
    }
  }
  massprofile=sumtot/(G*solmass)
  
  #return the caustic r200
  avg_density=massprofile/(4/3*pi*(ri[1:length(f_beta)])^3)
  finterp=approxfun(avg_density, ri[1:length(f_beta)])
  r200_est=finterp(200*crit)
  r500_est=finterp(500*crit)
  
  if(200*crit > max(avg_density, na.rm = T)){ # NEW
    warning('200 times critical density can not be interpolated', immediate. = T, call. = F)
    r200_est=r200
  }
    
  M200=tail(massprofile[ri[1:length(f_beta)] <= r200], 1)
  lout=list(r200_est=r200_est, r500_est=r500_est, M200=M200)
  return(lout)
}

run_caustic=function(gal_r=NA, gal_v=NA, clus_z=NA, r200=NA, clus_vdisp=NA, gal_memberflag=NA, rlimit=4, 
                     vlimit=3500, q=10, xmax=6, ymax=5000, cut_sample=T, mirror=T, H0=100, Om=.25, Ol=.75,
                     edge_int_remove=F, edge_perc=0.1, fbr=0.65, gapper=T, plot=T, verbose=T){
  
  r=gal_r; v=gal_v
  # calculate H(z)
  Hz = H0*sqrt(Om*(1+clus_z)^3 + Ol)
  hz=Hz/100 #little h(z)
  
  if(is.na(gal_memberflag)){
    data_set=data.frame(r=r, v=v)
  }else{
    data_set=data.frame(r=r, v=v, member=memberflags)
  }
  data_table=data_set
  
  #reduce sample within limits
  if(cut_sample) data_set=subset(data_set, r < rlimit & abs(v) < vlimit)
  
  if(gapper){
    data_set_copy=data_set
    class=shiftgapper(data_set[,1], data_set[,2], npbin = 25, plot = F)
    data_set=data_set[class == 0,]
    
    if(length(which(data_set[,1] < 2)) < 10){
      warning('Shiftgapper removing too much elements, skipping', immediate. = T, call. = F)
      data_set=data_set_copy
    }
  }
  
  if(nrow(data_set) < 2) stop('Data set has too few elements. Check the r and v objects. 
                              Could indicate wrong cluster/galaxy positions or redshifts')
  
  if(verbose) message(paste('DATA SET SIZE ', nrow(data_set)))
  
  rproj=data_set[,1]; vproj=data_set[,2]
  
  if(is.na(r200)){
    vdisp_prelim=biweightScale(vproj[rproj < 3], 9)
    r200_mean_prelim = 0.002*vdisp_prelim + 0.40
    r200 = r200_mean_prelim/1.7
    if(r200 > 3) r200 = 3.0
    rlimit=ifelse(r200 < 2, 3*r200, 5.5)
  }else{
    if(r200 > 3) r200 = 3.0
  }
  rlimit=max(rproj)
  #rlimit=ifelse(r200 < 2 & max(data_set[,1]) < 2*r200, 2*r200, max(data_set[,1])+.5)
  if(verbose) message(paste('Pre_r200 =', round(r200, 5)))
  
  # calculating density of phase-space -------------------------------------------------------------------
  if(mirror){
    if(verbose) message('Calculating Density w/Mirrored Data')
    gk=gaussian_kernel(rep(rproj,2), c(vproj,-vproj), r200, Hz, q, xmax, ymax, plot = F)
  }else{
    if(verbose) message('Calculating Density')
    gk=gaussian_kernel(rproj, vproj, r200, Hz, q, xmax, ymax, plot = F)
  }
  x_range=gk$x
  y_range=gk$y
  Zi=gk$z
  img_tot=Zi/max(abs(Zi))
  
  if(is.na(clus_vdisp)){
    v_cut=vproj[which(rproj < r200 & abs(vproj) < vlimit)]
    if(length(v_cut) < 2) v_cut=vproj[which(rproj < 2.5*r200 & abs(vproj) < vlimit)]
    pre_vdisp2=biweightScale(v_cut, 9)
    if(verbose) message(paste('Vdisp from galaxies = ', round(pre_vdisp2, 5)))
    pre_vdisp_comb=pre_vdisp2
  }else{
    pre_vdisp_comb=clus_vdisp
  }
  #if(verbose) message(paste('Combined Vdisp = ', round(pre_vdisp_comb, 5)))
  
  # Identify initial caustic surface and members within the surface --------------------------------------
  if(verbose) message('Calculating initial surface')
  if(is.na(gal_memberflag)){
    surface=findsurface(data_set, x_range, y_range, img_tot, r200, halo_vdisp=pre_vdisp_comb, mirror=mirror, 
                        edge_perc=edge_perc, Hz=Hz, edge_int_remove=edge_int_remove, q=q, plot = plot,
                        verbose = verbose, rimax=rlimit)
  }else{
    surface=findsurface(data_set, x_range, y_range, img_tot, r200, halo_vdisp=pre_vdisp_comb, mirror=mirror, 
                        edge_perc=edge_perc, Hz=Hz, edge_int_remove=edge_int_remove, q=q, plot = plot,
                        verbose = verbose, memberflags = data_set[,3], rimax=rlimit)
  }
  caustic_profile=surface$caustic_profile
  caustic_fit=surface$caustic_fit
  gal_vdisp=surface$gal_vdisp
  memflag=surface$memflag
  
  # caustic mass determination ---------------------------------------------------------------------------
  # Estimate the mass based off the caustic profile, beta profile (if given), and concentration (if given)
  if(!is.na(clus_z)){
    Mass=MassCalc(x_range, caustic_profile, gal_vdisp, clus_z, r200, fbr=fbr, H0=H0, Om=Om, Ol=Ol)
    r200_est=Mass$r200_est
    M200_est=Mass$M200
    if(verbose) message(paste('r200 estimate: ', round(r200_est, 5), 'Mpc'))
    if(verbose) message(paste('M200 estimate: ', format(M200_est, digits = 7), 'Msol'))
    Ngal=length(rproj[memflag == 1 & rproj <= r200_est])
  }
  
  # calculate velocity dispersion ------------------------------------------------------------------------
  vdisp_gal=biweightScale(vproj[memflag == 1])
  if(verbose) message('Vdisp estimate: ', round(vdisp_gal, 5), ' km/s')
  
  # defining membership ----------------------------------------------------------------------------------
  ids=as.numeric(rownames(data_set[memflag == 1,]))
  fit_outliers=caustic_outliers=rep(1, length(gal_r)); fit_outliers[ids]=0
  cinterp=approxfun(x_range, caustic_profile)
  vint=cinterp(rproj)
  mbr=ifelse(abs(vproj) <= vint, T, F)
  ids2=as.numeric(rownames(data_set[mbr,]))
  caustic_outliers[ids2]=0
  
  lout=list(data_table=data_table, data_set=data_set, x_range=x_range, y_range=y_range, img_tot=img_tot, 
            caustic_profile=caustic_profile, caustic_fit=caustic_fit, gal_vdisp=gal_vdisp,
            memflag=memflag, Ngal=Ngal, vdisp_gal=vdisp_gal, r200_est=r200_est, M200_est=M200_est,
            fit_outliers=fit_outliers, caustic_outliers=caustic_outliers)
  return(lout)
}

####################################### DEBUGGING ####################################### 
# r=a$dproj; v=a$vlos, clus_z=median(a$z)  # a is the cluster data
# r200=NA; clus_vdisp=NA; gal_memberflag=NA; rlimit=4; vlimit=3500; q=10; H0=100; xmax=6; ymax=5000 
# cut_sample=T; mirror=T; edge_int_remove=F; edge_perc=0.1; fbr=0.65; gapper=T; plot=T; verbose=T

### gaussian_kernel
# dproj=rep(data_set[,1],2); vlos=c(data_set[,2],-data_set[,2])
# norm=100; scale=10; xmax=6; ymax=5e3; by=0.05

### findsurface
# data=data_set; ri=x_range; vi=y_range; Zi=img_tot; halo_vdisp=pre_vdisp_comb; rimax=rlimit
# halo_scale_radius=NA; halo_scale_radius_e=0.01; bin=NA; beta=NA; maxv=5e3; memberflags=NA

### MassCalc
# ri=x_range; A=caustic_profile; vdisp=gal_vdisp; conc1=NA
