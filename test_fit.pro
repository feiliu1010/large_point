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

function myfunct,X,P
	amp=[0.1,1,10,0.1,0.1]
        ;P=[x0,E,B,Xp,theta]
        ex=exp(-(X-P[3]/amp[3])/(P[0]/amp[0]))
        ;define an array to truncate ex
        zero=fltarr(N_elements(X))
        zero[where(X ge P[3]/amp[3])]=1
        ex=ex*zero
        model= P[1]/amp[1]*GAUS_CONVOL(X,ex,P[4]/amp[4],nsigma=1.5)+P[2]/amp[2]
	return,model
end

pro test_fit
	amp=[0.1,1,10,0.1,0.1]
;	x_axle=[283.267,242.122,256.277,215.234,229.340, 188.320,202.413,161.495,$
;		175.540,134.644,148.677,107.882,121.867,81.0952,95.0673,54.3966,$
;		68.3214,27.6725,41.5845,1.03704,-12.8384,-25.6240,-39.4526,-52.1964,$
;		-66.0127,-78.7943,-92.5641,-105.304,-119.061,-131.838,-145.550,-158.285,$
;		-171.984,-184.757,-198.410,-211.140,-224.782, -237.549,-251.145,-263.869]
	x_axle=[-100,-85,-65,-50,-35,-10,0,20,40,55,65,85,105,125,145,160,175,195]
	x_axle1=x_axle[sort(x_axle)]
	upwind=-100
        downwind=200
        x_axle1=x_axle1[where( (x_axle1 le downwind) and (x_axle1 gt upwind) )]
;	Pin=[2.1310086/amp[0],11.254405/amp[1],12.996175/amp[2],2.7603894/amp[3],4.8022090/amp[4]]
;	Pin=[57.6,16.2D*10000,0.27*100000,15,21]
	Pin=[57.6,187D*6.02/(4.5*100),0.27,15,21]
	print,'pin',pin
	ex1=exp(-(x_axle1-Pin[3])/Pin[0])
;       define an array to truncate ex
        zero1=fltarr(N_elements(x_axle1))
        zero1[where(x_axle1 ge Pin[3])]=1
        ex1=ex1*zero1
	LD1=Pin[1]*GAUS_CONVOL(x_axle1,ex1,Pin[4],nsigma=3)+Pin[2]
	PLOT,x_axle1,LD1

	x_axle2=x_axle1
	LD2=[1.15886,1.10452,1.11014,1.21654,1.16865,1.41959,1.35578,1.68173,1.50582,$
	     1.99994,1.81950,2.53405,2.15938,2.92974,2.65980,3.89184,3.35491,4.67447,$
	     4.45565,4.18799,3.46221,2.49809,2.08695,1.69916,1.52623,1.49366,1.32464,$
	     1.27811,1.28128,1.30987,1.25986,1.18889,1.23028,1.30743,1.16301,1.16112,$
	     1.07913,1.00050,0.984788,1.04815]
	LD2=LD2[sort(x_axle)]
	distance=x_axle2[where(LD2 eq max(LD2))]
	print,'distance',distance

        ERR1=make_array(N_elements(x_axle1),value=0.1)
        weights1 = make_array(N_elements(x_axle1),value=1.0)
	ERR2=make_array(N_elements(x_axle2),value=0.1)
        weights2 = make_array(N_elements(x_axle2),value=1.0)
        ;Provide an initial guess of the function's parameters.
        ;start = [15.D,5.D,0.1D,5.1D,50.D]
	amp=[0.1,1,10,0.1,0.1]
        parinfo = replicate({value:0.D, fixed:0, limited:[0,0], $
                       limits:[0.D,0],relstep:0.D,mpside:0}, 5)
	parinfo[1].limited[0] = 1
        parinfo[1].limits[0]  = 0.D*amp[1]
        parinfo[0].limited[*] = 1
        parinfo[0].limits[0]  =10.D*amp[0]
        parinfo[0].limits[1]  =100.D*amp[0]
        parinfo[2].limited[*] = 1
        parinfo[2].limits[0]  =0.D*amp[2]
        parinfo[2].limits[1]  =(min(LD2)+10.D)*amp[2]
        parinfo[3].limited[*] = 1
        parinfo[3].limits[0]  =(distance-10.D)*amp[3]
        parinfo[3].limits[1]  =(distance+10.D)*amp[3]
        parinfo[4].limited[0] = 1
        parinfo[4].limits[0]  = 0.D*amp[4]
        
	parinfo[*].RELSTEP=0.05
	parinfo[*].mpside=2
	parinfo[*].value = [50.D*amp[0],max(LD2)*amp[1],min(LD2)*amp[2],distance*amp[3],20.D*amp[4]]
        print,'start',parinfo[*].value
	P1 = MPFITFUN('myfunct',x_axle1,LD1,ERR1,start,WEIGHTS=weights1,PARINFO=parinfo,/quiet)
	P2 = MPFITFUN('myfunct',x_axle2,LD2,ERR2,start,WEIGHTS=weights2,PARINFO=parinfo,yfit=yfit,/quiet)
;	print,'P1:',P1
	print,'P2:',P2

;	x_axle1_new=INTERPOLATE(x_axle1,findgen(100)*0.2)
	x_axle1_new=x_axle1
        F1=fltarr(N_elements(x_axle1_new))
        ex=exp(-(x_axle1_new-P1[3]/amp[3])/(P1[0]/amp[0]))
        ;define an array to truncate ex
	zero=fltarr(N_elements(x_axle1_new))
        zero[where(x_axle1_new ge P1[3]/amp[3])]=1
        ex=ex*zero
;	print,'ex1',ex
;	print,'Gx1',GAUS_CONVOL(x_axle1_new,ex,P1[4]/amp[4],nsigma=1.4)
	F1=P1[1]/amp[1]*GAUS_CONVOL(x_axle1_new,ex,P1[4]/amp[4],nsigma=1.4)+P1[2]/amp[2]

;	x_axle2_new=INTERPOLATE(x_axle2,findgen(N_elements(x_axle2))*1)
	x_axle2_new= x_axle2
	F2=yfit
        ;P=[x0,E,B,Xp,theta]
	x_axle=x_axle2
	x_axle_new=x_axle2_new
	LD=LD2
	F=F2

	PLOT,x_axle,LD,psym=2,$
        title='east & west',xtitle='x(km)',ytitle='NO2 Line Density(10^23 molec/cm)'
	OPLOT,x_axle_new,F
;	OPLOT,x_axle1,LD1,psym=4
;	OPLOT,x_axle1,F1
	image = tvrd(true =1)
	write_jpeg,'test.jpg',image,true=1

end
