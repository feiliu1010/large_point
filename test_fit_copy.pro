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
        model= P[1]/amp[1]*GAUS_CONVOL(X,ex,P[4]/amp[4],nsigma=1.4)+P[2]/amp[2]
	return,model
end

pro test_fit
	x_axle1=(findgen(20)-7)*15
	Pin=[4.4*3.1*3.6,187/(24*6.02),0.27D,15.D,21D]
	ex1=exp(-(x_axle1-Pin[3])/Pin[0])
;       define an array to truncate ex
        zero1=fltarr(N_elements(x_axle1))
        zero1[where(x_axle1 ge Pin[3])]=1
        ex1=ex1*zero1
	LD1=Pin[1]*GAUS_CONVOL(x_axle1,ex1,Pin[4],nsigma=2.5)+Pin[2]
;	x_axle=(findgen(30)-15)*13
	x_axle2=-[ -195.80911,      -177.77387,      -159.68170,      -141.64555,      -123.55256,$
      -105.51568,      -87.422055,      -69.384649,      -51.290599,      -33.252864,$
      -15.158582,       2.8792847,       20.973600,       39.011402,       57.105553,$
       75.143093,       93.236882,       111.27396,       129.36719,       147.40362,$
       165.49609,       183.53166]
	LD2=[0.114264,     0.109934,     0.118977,     0.149750,     0.144518,     0.172254,$
     0.221318,     0.312881,     0.399781,     0.472757,    0.523744,     0.454049,     0.380693,$
     0.344972,     0.292596,     0.253509,     0.181766,     0.135307,     0.107213,    0.0823060,    0.0780757,    0.0788648]*10
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
        parinfo[2].limits[1]  =2.D*amp[2]
        parinfo[3].limited[*] = 1
        parinfo[3].limits[0]  =0.D*amp[3]
	parinfo[3].limits[1]  =30.D*amp[3]
        parinfo[4].limited[0] = 1
        parinfo[4].limits[0]  = 0.D*amp[4]

	parinfo[*].RELSTEP=0.05
	parinfo[*].mpside=2
	parinfo[*].value = [50.D*amp[0],max(LD2)*amp[1],min([min(LD2),30])*amp[2],distance*amp[3],20.D*amp[4]]
        print,'start',parinfo[*].value
	P1 = MPFITFUN('myfunct',x_axle1,LD1,ERR1,start,WEIGHTS=weights1,PARINFO=parinfo,/quiet)
	P2 = MPFITFUN('myfunct',x_axle2,LD2,ERR2,start,WEIGHTS=weights2,PARINFO=parinfo,/quiet)
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
	F2=fltarr(N_elements(x_axle2_new))
        ;P=[x0,E,B,Xp,theta]
        ex=exp(-(x_axle2_new-P2[3]/amp[3])/(P2[0]/amp[0]))
        ;define an array to truncate exI
        zero=fltarr(N_elements(x_axle2_new))
        zero[where(x_axle2_new ge P2[3]/amp[3])]=1
        ex=ex*zero
        F2=P2[1]/amp[1]*GAUS_CONVOL(x_axle2_new,ex,P2[4]/amp[4],nsigma=1.4)+P2[2]/amp[2]

	x_axle=x_axle2
	x_axle_new=x_axle2_new
	LD=LD2
	F=F2

	PLOT,x_axle,LD,psym=2,$
        title='east & west',xtitle='x(km)',ytitle='NO2 Line Density(10^23 molec/cm)'
	OPLOT,x_axle_new,F
	OPLOT,x_axle1,LD1,psym=4
	OPLOT,x_axle1,F1
	image = tvrd(true =1)
	write_jpeg,'test.jpg',image,true=1

end
