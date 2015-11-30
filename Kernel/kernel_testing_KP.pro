
function gaussian_fit_1d,x,p

  y=p[0]*exp(-(x^2.0)/(2.0*((p[1]/SQRT(8.0*ALOG(2.0)))^2.0)))

  return,y

end


pro kernel_testing

mconv=mrdfits('~/san/Kernels/model_450c850_MK.fits',0,head)
bconv=mrdfits('~/san/Kernels/model_450c850_BK.fits',0,head)


model1=psf_gaussian(npixel=501,fwhm=13.0,/normalize)
model2=psf_gaussian(npixel=501,fwhm=48.0,/normalize)
model=(0.75*model1)+(0.25*model2)

omod1=psf_gaussian(npixel=501,fwhm=7.9,/normalize)
omod2=psf_gaussian(npixel=501,fwhm=25.0,/normalize)
omodel=(0.6*omod1)+(0.4*omod2)

x=findgen(501)-250

set_plot,'ps'
device,filename='~/san/Kernels/model_scuba2_beam_testing.ps',$
       xsize=20,ysize=15
!p.thick=3
!x.thick=3
!y.thick=3
!p.charthick=3


cgplot,x,model[250,*],color='blue',/ylog,xrange=[-50,50],$
       xtitle='Arcseconds',ytitle='Arbitrary Units',charsize=1.0

cgplot,x,omodel[250,*],color='dark grey',/overplot

cgplot,x,mconv[250,*],color='black',/overplot
cgplot,x,bconv[250,*],color='red',/overplot

cgplot,x,model[250,*],color='blue',/overplot

start=[0.003d,10.0d]

mfit=mpfitfun('gaussian_fit_1d',x[240:260],mconv[250,240:260],0.01*mconv[250,240:260],$
              start,perror=merrors,dof=mfree,bestnorm=mchi)

mmod=gaussian_fit_1d(x,mfit)
cgplot,x,mmod,/overplot,color='green',linestyle=5

bfit=mpfitfun('gaussian_fit_1d',x[240:260],bconv[250,240:260],0.01*mconv[250,240:260],$
              start,perror=merrors,dof=mfree,bestnorm=mchi)

bmod=gaussian_fit_1d(x,bfit)
cgplot,x,bmod,/overplot,color='orange',linestyle=5

cgtext,0.15,0.85,'Model 850 beam: 13"+48", eff. 14.1"',/normal,color='blue',charsize=1.0
cgtext,0.15,0.80,'Model 450 beam: 7.9"+25", eff. 9.6"',/normal,color='dark grey',charsize=1.0
cgtext,0.15,0.75,'Model-Kernel 850 beam: '+string(mfit[1],'(F0.2)')+'"',/normal,color='black',charsize=1.0
cgtext,0.15,0.70,'Beam-Kernel 850 beam: '+string(bfit[1],'(F0.2)')+'"',/normal,color='red',charsize=1.0


device,/close














END
