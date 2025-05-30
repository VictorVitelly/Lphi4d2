set terminal qt size 1900,1000
set lmargin 24
#set xlabel 'ϕ' font ',20'
set ylabel 'Frec.' font ',24' rotate by 0 offset -5,0
#set title 'Histogram for λ=0, μ_0^2=1' font ',20'
set key font ',22'
set xtics font ',17'
set ytics font ',17'

binwidth=3./101.
bin(x,width)=width*floor(x/width)
set boxwidth binwidth

a=1
sigma=1
xmean=0.
gauss(x)=a*exp(-(x-mean)**2 /(2*k**2)  )
fit gauss(x) '../data/histogram64.dat' u 1:2:3 via a, k, mean
chi2=(FIT_STDFIT*FIT_STDFIT)

set multiplot layout 4,1 title "Histogramas" font ',24'

        set yrange [0:0.95]

        plot '../data/histogram64.dat' u 1:2:3 notitle w errorbars, gauss(x) lt rgb "red" lw 2 title sprintf("Ajuste gaussiano: σ=%.4f±%.4f, <ϕ>=%.5f±%.5f, χ^2/dof=%.2f",k,k_err,mean,mean_err,chi2), '../data/histogram64.dat' u 1:2:3 smooth freq with boxes lt rgb "purple" lw 2 title "λ=0, μ_0^2=1"

        unset yrange
        binwidth=3.6/101.
        bin(x,width)=width*floor(x/width)
        set boxwidth binwidth

        plot '../data/histogram64.dat' u 4:5:6 notitle w errorbars, '../data/histogram64.dat' u 4:5:6 smooth freq with boxes lt rgb "purple" lw 2 title "λ=1, μ_0^2=-1.0"

        plot '../data/histogram64.dat' u 7:8:9 notitle w errorbars, '../data/histogram64.dat' u 7:8:9 smooth freq with boxes lt rgb "purple" lw 2 title "λ=1, μ_0^2=-1.4"

        binwidth=4.6/101.
        bin(x,width)=width*floor(x/width)
        set boxwidth binwidth
        set xlabel 'ϕ' font ',24'

        plot '../data/histogram64.dat' u 10:11:12 notitle w errorbars, '../data/histogram64.dat' u 10:11:12 smooth freq with boxes lt rgb "purple" lw 2 title "λ=1, μ_0^2=-2.0"

pause -1

unset multiplot
