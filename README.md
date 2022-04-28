# SPtools

SPtools merupakan kumpulan tools yang dibuat untuk mendukung pendugaan stok yang termasuk didalam kelompok surplus production. Tool pemodelan yang tersedia disini berfungsi:

1. Surplus produksi dengan asumsi equilibrium

Tool ini menghitung jumlah tangkapan ikan lestari (MSY) dan upaya penangkapan ikan lestari (Emsy) dengan asumsi equilibrium untuk model Schaefer dan Fox. Pendekatan ini ditampilkan disini hanya untuk tujuan edukasi sebagai contoh model yang akan memberikan estimasi MSY dan Emsy yang lebih tinggi, sehingga sangat tidak disarankan untuk dijadikan sebagai panduan dalam pengambilan kebijakan perikanan. Overestimasi reference point pada kondisi ketika status perikanan sedang dalam kondisi overexploited akan memberikan ilusi bahwa stok ikan masih banyak, sehingga dapat merugikan pelaku perikanan karena jumlah tangkapan yang rendah dan merugikan stok ikan karena semakin tingginya pemanfaatan. Banyak dari kita yang masih menggunakan metode ini, meskipun sudah tidak disarankan untuk digunakan sejak 1980an.

2. Surplus produksi dengan asumsi non equilibrium menggunakan metode time series fitting

Tool ini melakukan estimasi parameter K, B0, r, q dan menentukan jumlah tangkapan ikan lestari (MSY), biomassa ikan lestari (Bmsy), serta upaya penangkapan ikan lestari (Emsy) menggunakan data runut waktu dengan asumsi non-equilibrium untuk model Schaefer dan Fox. Tool ini sudah disesuaikan untuk kebutuhan data yang terbatas (dapat mengakomodasi ketiadaaan input data upaya penangkapan) serta sudah memperhitungkan kesalahan dalam pengambilan data (observation error) dan kesalahan dalam model (model error). Metode time series fitting disebut sebagai metode yang lebih baik dibandingkan dengan dua metode lain (metode equilibrium dan regresi) yang digunakan untuk melakukan estimasi parameter dalam model surplus produksi. Sebagian kecil dari kita sudah menggunakan metode ini, tetapi masih kurang tepat dalam melakukan analisisnya sehingga berakibat pada kurang tepatnya perhitungan MSY, Bmsy dan Emsy

3. Menghitung rentang data reference point

Tool ini melakukan estimasi atas reference point untuk jumlah tangkapan ikan lestari (MSY) dan upaya penangkapan ikan lestari (Emsy) serta menghitung standard error menggunakan data runut waktu dengan asumsi non-equilibrium untuk model Schaefer dan Fox. Pendekatan ini akan memberikan rentang kemungkinan error atas pendugaan MSY dan Emsy yang dilakukan.

4. Menyediakan data prior untuk parameter pertumbuhan r

Data prior setiap spesies untuk parameter pertumbuhan r diambil dari database fishbase dan sealifebase untuk mendukung pendugaan stok di tingkat spesies

5. Menghitung reference point untuk pengelolaan di tingkat spesies

Metode ini digunakan untuk menentukan MSY, Bmsy dan Emsy menggunakan pendekatan bayesian dan data runut waktu dengan asumsi non-equilibrium untuk model Schaefer dan Fox sehingga menghasilkan pendugaan stok yang lebih akurat di tingkat spesies

6. Membuat proyeksi atas kebijakan reference point berdasar tingkat pemanfaatan perikanan

Tool ini akan membuat grafik proyeksi biomass per biomass at msy (B/Bmsy) dan fishing per fishing at msy (F/Fmsy) sebagai panduan untuk melihat kebijakan yang dibuat saat ini serta memperkirakan limit dan target reference point


### Cara menggunakan SPtools untuk Surplus produksi dengan asumsi non equilibrium menggunakan data time series ####
#### membuat input data
df <- data.frame(tahun=c(...),
                 tangkapan=c(...),
                 upaya=c(...))

# SPtools sedang dalam proses penyempurnaan. Kami akan segera tanpilkan function yang ada dalam tool ini dalam sebuah package yang akan lebih mempermudah penggunaan model surplus production ini
