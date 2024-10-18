library(roahd)
library(fields)
library(viridis)
library(geoR)
library(laGP)
library(fda.usc)
library(geofd)
library(Rlinsolve)
library(progress)

###--- this function generate the spatial grid
##-- long   := a couple of min and max for the longitude
##-- lat    := a couple of min and max for the longitude 
##-- n_long := number of point for longitude
##-- n_lat  := number of point for latitude
##-- return a data frame with the coordinates of point
gene_grid = function(long,
                     lat = long,
                     n_long,
                     n_lat = n_long)
{
  n_long = n_long - 1
  n_lat = n_lat - 1
  # Defining the grid (spatial configuration of points)
  xlim_min = long[1]
  xlim_max = long[2]
  xinc = abs(long[2] - long[1]) / n_long
  ylim_min = lat[1]
  ylim_max = lat[2]
  yinc = abs(lat[2] - lat[1]) / n_lat
  par(mar = c(5, 5, 4, 1))
  grid <-
    expand.grid(seq(xlim_min, xlim_max, by = xinc),
                seq(ylim_min, ylim_max, by = yinc))
  grid = data.frame("Latitudine" = grid[, 1], "Longitudine" = grid[, 2])
  #plot(grid, xlab="x coordinate", ylab="y coordinate", col=1, pch=16, cex.axis=1.3,cex.lab=1.5)
  return("grid" = grid)
}
###--- this function compute the spatio-covariance matrix
##-- h    := the spazial distance between points
##-- ni   := residual variance
##-- c    := rate of decrease of spatial covariance
##-- return a data frame with the coordinates of point
covS = function(h, ni = 0.1, c = 0.1) {
  (1 - ni) * exp(-c * abs(h)) + ni * diag(nrow(h))
}

###--- this function generate the functional data
##-- G        := the functional mean
##-- E        := the functional error
##-- grid     := spatial coordinates
##-- scenario := the name of scenario
##-- return the functional data
simulated_data <- function(G, E, grid, scenario) {
  switch(
    scenario,
    "stat_sim" = G + E,
    "cubic_sim" = G ^ 3 + E,
    "skew_sim" = {
      X <- qgamma(pnorm(G, sd = sqrt(3)), 1, 1 / sqrt(3))
      X + E
    },
    "prod_sim" = sqrt(3) * G * abs(E),
    "gfun_sim" = {
      g <- function(x, p)
        sign(x) * abs(x) ^ p
      X <- g(G, grid + 1)
      X + E
    },
    "eastwest_sim" = {
      w <- pnorm(grid, 0.5, 0.1)
      sqrt(w) * G / sqrt(3) + sqrt(1 - w) * E
    },
    "nugget_sim" = G + grid * E,
    "spike_sim" = {
      dist <- sqrt((0.5 - grid) ^ 2)
      G + 10 * exp(-50 * (dist ^ 2))
    }
  )
}

###--- functional ordinary kriging
##--new.coords    := coordinates of new site
##-- coords       := coordinates of data
##-- data         := functional data
##-- smooth.type  := type of basis
##-- nbasis       := dimension of basis
##-- argvals      := interval
##-- lambda       := smoothing paramiter
##-- inputs for trace-variogram (see the help fit.tracevariog)
##-- L2norm       := distance of functional data with norm L2
##-- coefD        := coefficients of expanded basis
okfd2 = function (new.coords,
                  coords,
                  data,
                  smooth.type = NULL,
                  nbasis = max(50, dim(data)[1]),
                  argvals = seq(0, 1, len = dim(data)[1]),
                  lambda = 0,
                  cov.model = NULL,
                  fix.nugget = FALSE,
                  nugget = 0,
                  fix.kappa = TRUE,
                  kappa = 0.5,
                  max.dist.variogram = NULL,
                  L2norm,
                  coefD)
{
  smooth.type <- match.arg(smooth.type, c("bsplines",
                                          "fourier"))
  if (!is.null(cov.model)) {
    cov.model <- match.arg(cov.model,
                           c("spherical",
                             "exponential", "gaussian", "matern"))
  }
  if (is.null(new.coords))
    stop("new.coords is not an optional parameter")
  if (ncol(new.coords) != 2)
    stop("new.coords must be an n x 2 matrix")
  new.coords <- as.matrix(new.coords)
  coords <- as.matrix(coords)
  s <- dim(data)[2]
  
  # Empirical trace-variogram
  emp.trace.vari <- trace.variog(coords, L2norm)
  
  # Fitting a theoretical variogram model to empirical trace-variogram
  if (smooth.type == "fourier") {
    # partial sill is set to the quantile 0.75
    sigma2.0 <- quantile(emp.trace.vari$v, 0.75)
  } else{
    # partial sill is set to the variance of the data
    sigma2.0 <- var(as.vector(data))
  }
  trace.vari.objects <-
    fit.tracevariog(
      emp.trace.vari = emp.trace.vari,
      models = cov.model,
      sigma2.0 = sigma2.0,
      phi.0 = quantile(emp.trace.vari$Eu.d, 0.75),
      fix.nugget,
      nugget,
      fix.kappa,
      kappa,
      max.dist.variogram
    )
  
  fdmodel = list(emp.trace.vari = emp.trace.vari, trace.vari.objects = trace.vari.objects)
  a = fdmodel$trace.vari.objects$best
  #print(a$cov.model)
  trace.vari = fdmodel$trace.vari.objects$best
  Eu.d = fdmodel$emp.trace.vari$Eu.d
  s <- dim(coords)[1]
  
  # Distances among sites and distances to the prediction site
  # Distance matrix among sampling sites and distance to NEW SITES
  new.s <- dim(new.coords)[1]
  new.Eu.d <-
    as.matrix(dist(rbind(coords, new.coords), method = "euclidean"))
  new.Eu.d <-
    matrix(new.Eu.d[1:s, (s + 1):(s + new.s)], nrow = s, ncol = new.s)
  
  # Solving the system
  sigma2 <- trace.vari$cov.pars[1]
  leftmatrix <-
    sigma2 - cov.spatial(
      Eu.d,
      cov.model = trace.vari$cov.model,
      cov.pars = trace.vari$cov.pars,
      kappa = trace.vari$kappa
    )
  unosfila <- rep(1, s)
  leftmatrix <- rbind(leftmatrix, unosfila)
  unosycerocolumna <- c(rep(1, s), 0)
  leftmatrix <- cbind(leftmatrix, unosycerocolumna)
  
  rightmatrix <-
    sigma2 - cov.spatial(
      new.Eu.d,
      cov.model = trace.vari$cov.model,
      cov.pars = trace.vari$cov.pars,
      kappa = trace.vari$kappa
    )
  unosfila <- rep(1, new.s)
  rightmatrix <- rbind(rightmatrix, unosfila)
  
  functional.kriging.weights <-
    lsolve.cgs(leftmatrix, rightmatrix, verbose = F)[[1]]
  functional.kriging.weights.sinlagrange <-
    matrix(functional.kriging.weights[-(s + 1), ],
           nrow = s,
           ncol = new.s)
  
  # Solution
  
  krig.new.data <- coefD %*% functional.kriging.weights.sinlagrange
  
  # Prediction variance
  vect.semiv <- rightmatrix[-(s + 1), ]
  varianza <- functional.kriging.weights.sinlagrange * vect.semiv
  suma.varianza <- sum(varianza)
  pred.var <- suma.varianza + functional.kriging.weights[s + 1, ]
  
  prediction = list(
    pred = krig.new.data,
    var = pred.var,
    new.Eu.d = new.Eu.d,
    functional.kriging.weights = functional.kriging.weights
  )
  
  return.list <-
    list(
      coords = coords,
      data = data,
      argvals = argvals,
      nbasis = nbasis,
      lambda = lambda,
      new.coords = new.coords,
      emp.trace.vari = fdmodel$emp.trace.vari,
      trace.vari = fdmodel$trace.vari.objects$best,
      new.Eu.d = prediction$new.Eu.d,
      functional.kriging.weights = prediction$functional.kriging.weights,
      krig.new.data = prediction$pred,
      pred.var = prediction$var,
      trace.vari.array = fdmodel$trace.vari.objects$fitted,
      datafd = coefD
    )
  class(return.list) <- "geofd"
  return(return.list)
}

###--- conformal prediction
##-- Functions0   := functional data
##-- DOE0         := spatial coordinates
##-- argvals      := temporal grid
##-- fn           := number of band to compute
##-- alpha        := miss-coverage level
##-- dim_tr       := percentage of observations used as training set 
##-- modulazione  := type of modulation function used
##-- nncm         := type of non conformity measure used
##-- ns           := number of simulation
##-- return       := a list with predicted data, ray of prediction bans, moduletion function and quality index
convalida_cpk = function(Functions0,
                             DOE0,
                             argvals = seq(0, 1, length.out = nrow(Functions0)),
                             fn = 1,
                             alpha = 0.05,
                             dim_tr,
                             modulazione,
                             nncm,
                             ns = 0) {
  n <- dim(Functions0)[1]
  #  argvals <- Ti#seq(1,n, by=1)
  ######################################
  ### Cross validetion sulle basi
  ######################################
  #nb = 5:n
  #mdata = fdata(t(Functions0))
  #out<-optim.basis(mdata,lambda=0,numbasis=nb,type.basis="fourier")
  # out <- optim.basis(mdata,
  #                    lambda = 0,
  #                    numbasis = nb,
  #                    type.basis = "bspline")
  
  k = 30 #out$numbasis.opt
  zz = dim(DOE0)[1]
  b1.2 <-
    create.bspline.basis(rangeval = range(argvals) , nbasis = k)
  #b1.2 <- create.fourier.basis(rangeval =range(argvals), nbasis=K,period=1)#  fit the data without smoothing
  fd1.2 <- Data2fd(argvals = argvals, Functions0, basisobj = b1.2)
  plot(fd1.2)
  CfD0 = fd1.2$coefs
  M = getbasispenalty(b1.2, 0)
  L2ns0 = l2.norm(zz, fd1.2, M)
  bse = eval.basis(argvals, b1.2)
  lt = length(argvals)
  verifica = matrix(0, 1, zz) 
  verifica2l = matrix(0, 1, zz)
  raggi = matrix(0, 1, zz) #ray
  temp = matrix(0, nrow = zz , ncol = 1)
  Ps0s = matrix(0, lt, zz)
  As = matrix(0, k, zz)
  Ss = matrix(0, lt, zz)
  
  pb <- progress_bar$new(
    format = paste0("simulazione ", ns, " [:bar] :percent eta: :eta"),
    total = fn,
    clear = FALSE,
    width = 60
  )
  for (z in 1:fn) {
    s0 = z
    coord.cero <- matrix(DOE0[s0, ], nrow = 1, ncol = 2)
    rd = Functions0[, s0]
    
    Functions = Functions0[, -s0]
    DOE = as.matrix(DOE0[-s0, ])
    L2ns = L2ns0[-s0, -s0]
    CfD = CfD0[, -s0]
    
    n <- dim(Functions)[2]
    
    
    buffD = rbind(DOE0[s0, ], DOE)
    buffD = as.matrix(dist(buffD))
    buffD = buffD[1, -1]
    train = which(buffD < quantile(buffD, dim_tr))
    test = (1:n)[-train]
    l = length(test)
    lt = length(argvals)
    ## -- Conformal Prediction -  ##
    NCS = matrix(0, nrow = l, ncol = 1)
    
    pbuff = matrix(0, nrow = lt, ncol = l)
    start <- Sys.time()
    index_Crd = train
    index_X = index_Crd
    
    okfd.res <- okfd2(
      new.coords = coord.cero,
      coords = DOE[index_Crd, ],
      data = Functions[, index_X],
      nbasis = k,
      argvals = argvals,
      fix.nugget = TRUE,
      L2norm = L2ns[index_X, index_X] ,
      coefD = CfD[, index_X]
    )
    Ps0 = okfd.res$krig.new.data
    plot(okfd.res)
    
    for (j in (1:l)) {
      index_Crd = c(train, test[j])
      buff = Functions[, index_Crd]
      buff2 = CfD[, index_Crd]
      index_X = index_Crd
      okfd.res <- okfd2(
        new.coords = coord.cero,
        coords = DOE[index_X, ],
        data = buff,
        nbasis = k,
        argvals = argvals,
        fix.nugget = TRUE,
        L2norm = L2ns[index_X, index_X] ,
        coefD = buff2
      )
      pbuf = okfd.res$krig.new.data
      pbuff[, j] = bse %*% pbuf
      #      cat(" - ", j)
    }
    S = matrix(0, nrow = lt, ncol = 1)
    A = Ps0
    Ps0 = bse %*% Ps0
    switch(modulazione,
           "1" = for (j in 1:lt) {
             S[j] = max(abs((pbuff[j,] - Ps0[j])), na.rm = T)
           },
           "2" = {
             S = sqrt(rowMeans((
               pbuff - matrix(Ps0, nrow = nrow(pbuff), ncol = ncol(pbuff))
             ) ^ 2))
           })
    switch(nncm,
           "1" = for (j in 1:l) {
             NCS[j] = max(abs((pbuff[, j] - Ps0) / S), na.rm = T)
           },
           "2" = for (j in 1:l) {
             NCS[j] = mean((pbuff[, j] - Ps0) ^ 2 / S)
             NCS[j] = sqrt(NCS[j])
           })
    end <- Sys.time()
    #    cat("\n")
    #    print(end - start)
    if (1) {
      tempo_trasc = (end - start)
    } else{
      tempo_trasc = tempo_trasc + (end - start)
    }
    #    print(paste0("estimated end at ", Sys.time() + (zz - z) * tempo_trasc))
    temp[z] = tempo_trasc
    r = (quantile(NCS, 1 - alpha))
    low = Ps0 - r * S
    #rd=bse%*%CfD0[,z]
    a = low - rd
    up = Ps0 + r * S
    b = rd - up
    verifica[z] = ((length(a[a > 0]) != 0) ||
                     (length(b[b > 0]) != 0))
    raggi[z] = r
    As[, z] = A
    Ps0s[, z] = Ps0
    Ss[, z] = S
    verifica2l[z] = length(a[a > 0]) + length(b[b > 0])
    pb$tick()
  }
  
  
  ###
  ### verifica2lp è il cov alfa funzionale!
  ###
  covAG = mean(1 - verifica) * 100
  covAL = (1 - verifica2l / lt) * 100
  covAG
  covAL
  hist(covAL)
  
  
  for (i in 1:fn) {
    Ps0 = As[, i] %*% t(bse)
    r = raggi[i]
    S = Ss[, i]
    rd = Functions0[, i]
    low = Ps0 - r * S
    up = Ps0 + r * S
    plot(
      argvals,
      Ps0,
      col = 1,
      lwd = 2,
      ylim = range(c(up, low, rd)),
      type = "l",
      lty = 1,
      main = paste("Prediction - curve", i),
      xlab = "Day",
      ylab = "Temperature (Degrees C)"
    )
    lines(argvals, up, col = "red")
    lines(argvals, low, col = "blue")
    lines(
      argvals,
      rd ,
      type = "p",
      pch = 20,
      cex = 0.5,
      col = 3,
      lwd = 1
    )
    rd = bse %*% CfD0[, i]
    lines(
      argvals,
      rd ,
      type = "l",
      pch = 20,
      cex = 0.5,
      col = 3,
      lwd = 1
    )
  }
  return(list(argvals, As, raggi, Ss, covAG, covAL, temp, bse, b1.2))
}

###--- this function generate the scenario
##-- Mu       := the mean function
##-- err      := the error function
##-- nun_scen := the type of scenario
sceni = function(Mu, err, nun_scen) {
  # Plot dei dati funzionali simulati
  # Lista degli scenari
  scenarios <- c(
    "stat_sim",
    "cubic_sim",
    "skew_sim",
    "prod_sim",
    "gfun_sim",
    "eastwest_sim",
    "nugget_sim",
    "spike_sim"
  )
  
  # Ciclo per plottare tutti gli scenari
  #  par(mfrow = c(2, 4))  # Imposta la disposizione dei plot
  Scemy = list()
  for (scenario in scenarios[nun_scen]) {
    # Genera dati funzionali simulati per lo scenario corrente
    Data <- simulated_data(Mu, err, Ti, scenario)
    # Data = t(Data)
    # Plot dei dati funzionali simulati
    plot(
      fData(Ti, Data),
      main = paste("Simulazione -", scenario),
      xlab = "grid",
      lwd = 2
    )
    Scemy[[scenario]] = Data
  }
  par(mfrow = c(1, 1))
  return(Scemy)
}


