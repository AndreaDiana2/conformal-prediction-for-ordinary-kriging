source("function.R")
main_directory = getwd()
#
sub_directory_simu = "simulation_plan"
dir.create(sub_directory_simu)
setwd(sub_directory_simu)
##--- number of current simulation ---###
ns = 1

###--- dimension of spatial grid ---###
dim_sp_grig = matrix(c(10, 10),
                         ncol = 2,
                         nrow = 1,
                         byrow = T)

###--- number of observations ---###
dim_dati = c(100)

###--- percentage of observation for training set ---###
dimensione_train = c(0.25, 0.5, 0.75)

###--- parameter of covariance function, c---###
c = c(0.1, 0.9)
###--- parameter of covariance function, ni---###
ni = c(0.1, 0.9)

for (size_dati in 1:length(dim_dati)) {
  sub_directory_simu = getwd()
  c_data_dim = paste0("dim_dateset", dim_dati[size_dati])
  dir.create(c_data_dim)
  setwd(c_data_dim)
###--- create a subdirectory ---###  
  ###--- size_dati identifies the number of observation ---###  
  
  for (scenario in 1:2) {
    c_data_dim = getwd()
    c_scenario = paste0("scenario", scenario)
    dir.create(c_scenario)
    setwd(c_scenario)
###--- create a subdirectory ---###  
    ###--- scenario identifies the scenario ---###  
    
    for (simulazione in 1:1) {
      c_scenario = getwd()
      c_simul = paste0("simulazione", simulazione)
      dir.create(c_simul)
      setwd(c_simul)
###--- create a subdirectory ---###        
      ###--- simulazione identifies the number of simulation ---###        

      for (size_train in 1:3) {
        c_simul = getwd()
        c_train = paste0("dim_train", dimensione_train[size_train])
        dir.create(c_train)
        setwd(c_train)
###--- create a subdirectory ---###  
        ###--- size_train identifies the percentage of data to use as training set ---###  
        
        for (modulazione in 1:2) {
          c_train = getwd()
          c_modulazione = paste0("modulazione", modulazione)
          dir.create(c_modulazione)
          setwd(c_modulazione)
###--- create a subdirectory ---###  
          ###--- modulazione identifies the modulation function to use in the simulation  ---###  

          for (nncm in 1:2) {
            c_modulazione = getwd()
            c_nncm = paste0("nncm", nncm)
            dir.create(c_nncm)
            setwd(c_nncm)
###--- create a subdirectory ---###  
            ###--- nncm identifies the non conformity measure to use in the simulation  ---###  
            
            for (i1 in 1:2) {
###--- create a subdirectory ---###  
              ###--- i1 identifies the parameter c in the covariance function---###  
              
              for (i2 in 1:2) {
###--- create a subdirectory ---###  
                ###--- i2 identifies the parameter ni in the covariance function  ---###  
                
                set.seed(simulazione)
                nlo = dim_sp_grig[size_dati, 1]
                nla = dim_sp_grig[size_dati, 2]
                ###--- create the spatial grid 10x10 tra [-1,1]x[0,1]
                coord = gene_grid(
                  long = c(-1, 1),
                  lat = c(0, 1),
                  n_long = nlo,
                  n_lat = nla
                )
                plot(
                  coord,
                  xlab = "x coordinate",
                  ylab = "y coordinate",
                  col = 1,
                  pch = 16,
                  cex.axis = 1.3,
                  cex.lab = 1.5
                )
                ###--- spatial distance
                H = as.matrix(dist(coord))
                ###--- here we generate the functional data
                M = 100
                ###--- the temporal grid
                Ti = seq(0, 1, length.out = M) * pi * 2
                
                simf1 <- function(t) {
                  -2 * sin(t - 1) * log(t + 0.5)
                }
                ###--- covariance function
                CSH = covS(H, c[i1], ni[i2])
                N = dim(CSH)[1]
                Ti2 <-
                  seq(range(Ti)[1], range(Ti)[2], length.out = N)
                err =  t(generate_gauss_fdata(
                  M,
                  centerline = simf1(Ti2),
                  Cov = CSH
                ))
                m = 1 / 2 * (Ti / (2 * pi)) + 2 * sin(Ti)
                Mu = t(matrix(m, M, N))
                ###--- generate the scenario
                Sc = sceni(Mu, err, scenario)
                
                ## conformal prediction
                a = convalida_cpk(
                  t(Sc[[1]]),
                  coord,
                  Ti,
                  fn = N,
                  alpha = 0.05,
                  dimensione_train[size_train],
                  modulazione,
                  nncm,
                  ns
                )
                ns = ns + 1
                save.image(
                  paste0(
                    "res_nncm",
                    nncm,
                    "_mod",
                    modulazione,
                    "_dimTrain",
                    dimensione_train[size_train],
                    "_sim",
                    simulazione,
                    "c",
                    c[i1],
                    "ni",
                    ni[i2],
                    "_scen",
                    scenario,
                    "_dimData",
                    dim_dati[size_dati],
                    ".RData"
                  )
                )
              }
            }
            
            setwd(c_modulazione)
          }
          setwd(c_train)
        }
        setwd(c_simul)
      }
      setwd(c_scenario)
    }
    setwd(c_data_dim)
  }
  setwd(sub_directory_simu)
}

setwd(main_directory)