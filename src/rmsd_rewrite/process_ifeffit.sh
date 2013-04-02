head="
# read data\n
read_data(file=chi.chi3,type=chi,group=data)\n
my.k   = data.k\n
my.chi = data.chi/(data.k**3)\n
# Interpolate to grid\n
newdata.k = range(0,10.0,0.05)\n
newdata.chi = qinterp(my.k, my.chi, newdata.k)\n
# Guess e0 value\n
guess e  =  0.0\n
set   s  = 1.0\n
guess s1 = 0.0025\n
guess s2 = 0.0005\n
guess s3 = 0.0005\n
guess s4 = 0.0005\n
guess dr = 0.02\n
# Read Paths"

tail="
# define FFT Paramters\n
set (kmin = 1.0, kmax =10.0)\n
set (kweight=3,dk = 1, kwindow='hanning')\n
set (rmin = 0, rmax = 6)\n
# Create Initial Guess\n
ff2chi(1-100,group = init)\n
# Do the fit\n
feffit(1-100,chi = newdata.chi, group=fit )\n
fftf(newdata.chi)\n
# Write results\n
show @paths\n
show @variables\n
fit3.chi = fit.chi*fit.k**3\n
init3.chi = init.chi*init.k**3\n
write_data(file=my_chi.chi3,fit.k,fit3.chi)\n
write_data(file=fit_dw_init.DAT,init.k,init3.chi)\n
write_data(file=my_FT.out,fit.r,fit.chir_mag)\n
write_data(file=exp_r.out,newdata.r,newdata.chir_mag)\n
exit"

$(echo -e $head > process.iff)

COUNTER=1
for i in `ls feff00*`
do
	echo "path($COUNTER, file=$i, e0 = e, s02 = s, sigma2 = s1)" >> process.iff
	COUNTER=$[$COUNTER + 1]
done

$(echo -e $tail >> process.iff)