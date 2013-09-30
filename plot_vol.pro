function gaus_convol,x,y,sig,nsigma=nsigma

on_error,2
c=check_math(0,1)
if n_params() ne 3 then message,'usage: result=gaus_convol(x,y,sig)'

if n_elements(nsigma) ne 1 then nsigma=2.5      ; approximate width in
                                                ; units of sigma
nsigma2= NSIGMA*2
n      = n_elements(x)
conv   = (max(double(x))-min(double(x)))/(n -1) ; conversion, units/point
n_pts  = ceil(nsigma2*sig/conv)                 ; number of points
n_pts  = (n_pts > 2) < (n-2)                    ; restrictions on n_pts
if (n_pts/2*2) eq n_pts then n_pts = n_pts +1   ; make odd number of points
xvar = (dindgen(n_pts)/(n_pts-1)-0.5d0)*n_pts   ; approx. - NSIGMA < x < +NSIGMA

gaus=exp(-.5*(xvar/(sig/conv))^2)               ; gaussian of width sig.
return,convol(y,gaus,/center)/total(gaus)       ; do convolution.
end

pro plot_vol
        x_axle1=(findgen(20)-6)*13
        Pin=[4.4*3.1*3.6,187.D*6.02*3600.*3.1/10000000,0.27D,15.D,21D]
	print,'pin',pin
        ex1=exp(-(x_axle1-Pin[3])/Pin[0])
;       define an array to truncate ex
        zero1=fltarr(N_elements(x_axle1))
        zero1[where(x_axle1 ge Pin[3])]=1
        ex1=ex1*zero1
        LD1=Pin[1]*GAUS_CONVOL(x_axle1,ex1,Pin[4],nsigma=1)+Pin[2]

	x_axle1_new=INTERPOLATE(x_axle1,findgen(200)*0.1)
        Pin=[4.4*3.1*3.6,1.3D,0.27D,15D,21D]
        print,'pin2',pin
        ex1=exp(-(x_axle1_new-Pin[3])/Pin[0])
;       define an array to truncate ex
        zero1=fltarr(N_elements(x_axle1_new))
        zero1[where(x_axle1_new ge Pin[3])]=1
        ex1=ex1*zero1
        LD2=Pin[1]*GAUS_CONVOL(x_axle1_new,ex1,Pin[4],nsigma=1)+Pin[2]
	PLOT,x_axle1,LD1,psym=2,$
        title='east & west',$
	xtitle='x(km)',ytitle='NO2 Line Density(10^23 molec/cm)'
;	OPLOT,x_axle1_new,LD2
end
