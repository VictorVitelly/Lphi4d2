set lmargin 20
set xlabel 'ϕ' font ',20'
set ylabel 'Freq.' font ',20' rotate by 0
set title 'Histogram for λ=0, μ_0^2=1' font ',20'
set key font ',18'
binwidth=3./101.
bin(x,width)=width*floor(x/width)
set boxwidth binwidth
a=1
sigma=1
xmean=0.
gauss(x)=a*exp(-(x-mean)**2 /(2*k**2)  )
fit gauss(x) '../data/histogram.dat' u 1:3:5 via a, k, mean
chi2=(FIT_STDFIT*FIT_STDFIT)
#plot 'data/data.dat' using (bin($2,binwidth)):(1.0) smooth freq with boxes notitle, gauss(x)
plot '../data/histogram.dat' u 1:3:5 notitle w errorbars, gauss(x) lw 2 title sprintf("Gaussian fit: σ=%.4f±%.4f, <ϕ>=%.5f±%.5f, χ^2/dof=%.2f",k,k_err,mean,mean_err,chi2), '../data/histogram.dat' u 1:3:5 smooth freq with boxes notitle lt rgb "purple" lw 2
#plot 'data/data.dat' using (bin($2,binwidth)):(1.0) smooth freq with boxes notitle
pause -1
