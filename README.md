# SPtools

SPtools merupakan kumpulan tools yang dibuat untuk mendukung pendugaan stok yang termasuk didalam kelompok surplus production. Tool pendugaan stok yang tersedia disini adalah:

1) Surplus production model - Schaefer dan Fox dengan asumsi equilibrium
2) Surplus production model - Schaefer, non-equilibrium
3) Surplus production model - Fox, non-equilibrium
4) Rextractor to extract r prior from fishbase and sealifebase database
5) S.NEq_Projector
6) F.NEq_Projector

### Cara menggunakan package SPtools ####
# membuat input data
df <- data.frame(tahun=c(...),
                 tangkapan=c(...),
                 upaya=c(...))

# melakukan estimasi parameter menggunakan optim
K <- max(df[,2])*2
Bo <- max(df[,2])*1.5
r <- 0.5
q <- 0.00025

startPars <- c(log(K), log(Bo), log(r), log(q), log(0.1))

fit <- optim(par=startPars, 
             fn=Schaefer_NEq_NLL, 
             df=df, 
             method="Nelder-Mead")

# membuat plot hasil data fitting
predicted <- NEq.Sch(inpars = fitted_pars[1:4], df)

# membuat proyeksi
