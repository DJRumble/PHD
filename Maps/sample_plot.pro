pro sample_plot
;makes up a mock plot for Damian

set_plot,'ps'
device, filename='sample_plot.ps'
device, /color,bits_per_pixel=8
loadct, 13
!p.thick=4 & !x.thick=4 & !y.thick=4 & !p.charthick=4

;create two random distributions

random_x1 = randomu(seed,100)
random_x2 = randomn(seed,100)
random_y1 = randomu(seed,100)*0.5
random_y2 = randomn(seed,100)*0.5+0.5

;main plot
plot, random_x1, random_y1, xrange=[-1,2],yrange=[-1,2], xtitle='Random axis 1 (units)',$
	ytitle='Random axis 2 (units)', position=[0.2,0.2,0.7,0.7], /nodata
oplot, random_x1, random_y1, psym=4, color=50
oplot, random_x2, random_y2, psym=6, color=250

;basic histogram setup
binsz = 0.1
minbin = -1. & maxbin=2.
nbins = (maxbin-minbin)/binsz
bins = findgen(nbins)*binsz+minbin + binsz/2.


;Histograms for the x axis (put on top of main plot)
histx1 = histogram(random_x1,min=minbin,max=maxbin,binsize=binsz)
histx2 = histogram(random_x2,min=minbin,max=maxbin,binsize=binsz)

plot, bins, histx1, xtitle=' ', ytitle='Number', xrange=[-1,2], $
	yrange=[-0.2,max(histx1)+2.], xtickname=[' ',' ',' ',' ',' ',' ',' '],$
	position=[0.2,0.7,0.7,0.9],/nodata,/noerase
oplot, bins, histx1,psym=10, color=50
oplot, bins, histx2,psym=10, color=250


;Histograms for the y axis (put to the right of the main plot)
histy1 = histogram(random_y1,min=minbin,max=maxbin,binsize=binsz)
histy2 = histogram(random_y2,min=minbin,max=maxbin,binsize=binsz)

plot, histy1, bins, xtitle='Number', ytitle=' ', yrange=[-1,2], $
	xrange=[-0.5,max(histy1)+2.], ytickname=[' ',' ',' ',' ',' ',' ',' '],$
	position=[0.7,0.2,0.9,0.7],/nodata,/noerase
;more complicated for a sideways histogram
for i = 1, nbins-1 do begin
	;top of bin
	oplot, histy1(i)*[1,1],[bins(i-1),bins(i)], color=50
	oplot, histy2(i)*[1,1],[bins(i-1),bins(i)], color=250
	;left / top side of bin
	oplot, [histy1(i-1),histy1(i)],bins(i-1)*[1,1], color=50
	oplot, [histy2(i-1),histy2(i)],bins(i-1)*[1,1], color=250
endfor


device,/close
set_plot,'x'

stop

end
