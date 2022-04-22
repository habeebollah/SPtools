# SPtools

SPtools merupakan kumpulan tools yang dibuat untuk mendukung pendugaan stok yang termasuk didalam kelompok surplus production. Tool pemodelan yang tersedia disini berfungsi:

1. Surplus produksi dengan asumsi equilibrium
Tool ini menghitung jumlah tangkapan ikan lestari (MSY) dan upaya penangkapan ikan lestari (Emsy) dengan asumsi equilibrium untuk model Schaefer dan Fox. Pendekatan ini ditampilkan disini hanya untuk tujuan edukasi sebagai contoh model yang akan memberikan estimasi MSY dan Emsy yang lebih tinggi, sehingga sangat tidak disarankan untuk dijadikan sebagai panduan dalam pengambilan kebijakan perikanan.

2. Surplus produksi dengan asumsi non equilibrium menggunakan data time series
Tool ini melakukan estimasi parameter K, B0, r, q dan menentukan jumlah tangkapan ikan lestari (MSY), biomassa ikan lestari (Bmsy), serta upaya penangkapan ikan lestari (Emsy) menggunakan data runut waktu dengan asumsi non-equilibrium untuk model Schaefer dan Fox. Tool ini sudah disesuaikan untuk kebutuhan data yang terbatas (dapat mengakomodasi ketiadaaan input data upaya penangkapan) serta sudah memperhitungkan kesalahan dalam pengambilan data (observation error) dan kesalahan dalam model (model error).

3. Menghitung rentang data reference point
Tool ini melakukan estimasi atas reference point untuk jumlah tangkapan ikan lestari (MSY) dan upaya penangkapan ikan lestari (Emsy) serta menghitung standard error menggunakan data runut waktu dengan asumsi non-equilibrium untuk model Schaefer dan Fox. Pendekatan ini akan memberikan rentang kemungkinan error atas pendugaan MSY dan Emsy yang dilakukan.

4. Menyediakan data prior untuk parameter pertumbuhan r
Data prior setiap spesies untuk parameter pertumbuhan r diambil dari database fishbase dan sealifebase untuk mendukung pendugaan stok di tingkat spesies

5. Menghitung reference point untuk pengelolaan di tingkat spesies
Metode ini digunakan untuk menentukan MSY, Bmsy dan Emsy menggunakan pendekatan bayesian dan data runut waktu dengan asumsi non-equilibrium untuk model Schaefer dan Fox sehingga menghasilkan pendugaan stok yang lebih akurat di tingkat spesies

6. Membuat proyeksi atas kebijakan reference point berdasar tingkat pemanfaatan perikanan
Tool ini akan membuat grafik proyeksi biomass per biomass at msy (B/Bmsy) dan fishing per fishing at msy (F/Fmsy) sebagai panduan untuk melihat kebijakan yang dibuat saat ini serta memperkirakan limit dan target reference point


### Cara menggunakan package SPtools ####
#### membuat input data
df <- data.frame(tahun=c(...),
                 tangkapan=c(...),
                 upaya=c(...))

#### melakukan estimasi parameter menggunakan optim
K <- max(df[,2])*2
Bo <- max(df[,2])*1.5
r <- 0.5
q <- 0.00025

startPars <- c(log(K), log(Bo), log(r), log(q), log(0.1))

fit <- optim(par=startPars, 
             fn=Schaefer_NEq_NLL, 
             df=df, 
             method="Nelder-Mead")

#### membuat plot hasil data fitting
predicted <- NEq.Sch(inpars = fitted_pars[1:4], df)

#### membuat proyeksi
