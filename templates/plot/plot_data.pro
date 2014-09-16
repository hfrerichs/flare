;- grid ----------------------------------------------------------------
pro siscape, lun, sout, iout
	str	= ''
	readf, lun, str
	print, strmid(str,2,80)
	sout	= strmid(str,2,20)
	iout	= long(strmid(str,32,10))
end; siscape
pro iscape, lun, iout
	str	= ''
	readf, lun, str
	print, strmid(str,2,80)
	iout	= long(strmid(str,32,10))
end; iscape
pro rscape, lun, rout
	str	= ''
	readf, lun, str
	print, strmid(str,2,80)
	rout	= double(strmid(str,32,10))
end; rscape

pro read_grid, grid_file, xoffset=xoffset, yoffset=yoffset
	common	grid,	grid_id, x_data, y_data, n_x, n_y, n_xy, irr, isort, x_title, y_title, $
			x_min, x_max, y_min, y_max

	str	= ''
	openr, lun, grid_file, /get_lun
	readf, lun, str
	if (strmid(str,2,7) ne 'grid_id') then begin
		print, 'no grid type defined'
		stop
	endif
	grid_id	= fix(strmid(str,12,2))
	print, 'grid id = ', grid_id

	case grid_id of
	0: begin
		iscape, lun, n_xyz

		x_data	= make_array(n_xyz, /fl)
		y_data	= make_array(n_xyz, /fl)
		tmp	= make_array(3, /fl)
		for i=0L,n_xyz-1 do begin
			readf, lun, tmp
			x_data(i)	= tmp(0)
			y_data(i)	= tmp(2)
		endfor
		n_xy	= n_xyz
		isort	= 0
	end
	1: begin
		iscape, lun, n_RZ
		rscape, lun, phi

		x_data	= make_array(n_RZ, /fl)
		y_data	= make_array(n_RZ, /fl)
		x_title	= 'major radius [cm]'
		y_title	= 'Z [cm]'

		tmp	= make_array(2, /fl)
		for i=0L,n_RZ-1 do begin
			readf, lun, tmp
			x_data(i)	= tmp(0)
			y_data(i)	= tmp(1)
		endfor
		n_xy	= n_RZ
		isort	= 0
	end
	2: begin
		iscape, lun, n_rt
		rscape, lun, r_center
		rscape, lun, phi

		x_data	= make_array(n_rt, /fl)
		y_data	= make_array(n_rt, /fl)
		x_title	= 'poloidal angle [deg]'
		y_title	= 'minor radius [cm]'

		tmp	= make_array(2, /fl)
		for i=0L,n_rt-1 do begin
			readf, lun, tmp
			x_data(i)	= 360.0 - tmp(1)
			y_data(i)	= tmp(0)
		endfor
		n_xy	= n_rt
		isort	= 0
	end
	3: begin
		iscape, lun, n_R
		iscape, lun, n_Z

		x_data	= make_array(n_R, /fl)
		y_data	= make_array(n_Z, /fl)
		x_title	= 'R [cm]'
		y_title	= 'Z [cm]'

		rscape, lun, phi

		readf, lun, x_data
		readf, lun, y_data
		n_x	= n_R
		n_y	= n_Z
		n_xy	= n_x * n_y
		isort	= 1
	end
	4: begin
		iscape, lun, n_t
		iscape, lun, n_r

		x_data	= make_array(n_t, /fl)
		y_data	= make_array(n_r, /fl)
		x_title	= 'poloidal angle [deg]'
		y_title	= 'minor radius [cm]'

		rscape, lun, r_center
		rscape, lun, phi

		readf, lun, y_data
		readf, lun, x_data
		;x_data	= 360.0 - x_data
		n_x	= n_t
		n_y	= n_r
		n_xy	= n_x * n_y
		isort	= 1
	end
	41: begin;	same as 4 but with x and y coordinates switched for plotting
		iscape, lun, n_t
		iscape, lun, n_r

		x_data	= make_array(n_r, /fl)
		y_data	= make_array(n_t, /fl)
		x_title	= 'Minor Radius [cm]'
		y_title	= 'Poloidal Angle [deg]'

		rscape, lun, r_center
		rscape, lun, phi

		readf, lun, x_data
		readf, lun, y_data
		y_data	= 360.0 - y_data
		n_x	= n_r
		n_y	= n_t
		n_xy	= n_x * n_y
		isort	= -1
	end
	5: begin
		iscape, lun, n_p
		iscape, lun, n_t

		x_data	= make_array(n_p, /fl)
		y_data	= make_array(n_t, /fl)
		x_title	= 'toroidal angle [deg]'
		y_title	= 'poloidal angle [deg]'

		rscape, lun, r_center
		rscape, lun, r_min

		readf, lun, x_data
		readf, lun, y_data
		n_x	= n_p
		n_y	= n_t
		n_xy	= n_x * n_y
		isort	= 1
	end
	6: begin
		iscape, lun, n_RZ
		iscape, lun, n_p

		x_data	= make_array(n_p, /fl)
		y_data	= make_array(n_RZ, /fl)
		x_title	= 'toroidal angle [deg]'
		y_title	= 'length [cm]'

		R	= 0.0
		Z	= 0.0
		L	= 0.0
		for i=0,n_RZ-1 do begin
			readf, lun, R, Z, L
			y_data(i)	= L
		endfor
		readf, lun, x_data
		n_x	= n_p
		n_y	= n_RZ
		n_xy	= n_x * n_y
		isort	= 1
	end
	7: begin;	Theta-Psi grid (irregular)
		iscape, lun, n
		rscape, lun, phi

		rscape, lun, psi_axis
		rscape, lun, psi_sepx

		x_data	= make_array(n, /fl)
		y_data	= make_array(n, /fl)
		x_title	= 'poloidal angle [deg]'
		y_title	= 'normalized flux'

		theta	= 0.0
		psi	= 0.0
		for i=0L,n-1 do begin
			readf, lun, theta, psi
			x_data(i)	= theta
			y_data(i)	= psi
		endfor
		n_xy	= n
		isort	= 0
	end
	8: begin;	Theta-Psi grid (regular)
		iscape, lun, n_theta
		iscape, lun, n_psi
		rscape, lun, phi

		x_data	= make_array(n_theta, /fl)
		y_data	= make_array(n_psi, /fl)
		x_title	= 'poloidal angle [deg]'
		y_title	= 'normalized flux'

		rscape, lun, psi_axis
		rscape, lun, psi_sepx

		readf, lun, x_data
		readf, lun, y_data
		;x_data	= 360.0 - x_data
		n_x	= n_theta
		n_y	= n_psi
		n_xy	= n_x * n_y
		isort	= -1
	end
	; regular phase space grid
	50: begin
		iscape, lun, n_x
		iscape, lun, n_v

		x_data	= make_array(n_x, /fl)
		y_data	= make_array(n_v, /fl)
		x_title	= 'real space'
		y_title	= 'velocity space'

		readf, lun, x_data
		readf, lun, y_data
		n_y	= n_v
		n_xy	= n_x * n_y
		isort	= 1
	end
	; regular grid, user defined labels
	10: begin
		siscape, lun, x_title, n_x
		siscape, lun, y_title, n_y

		x_data	= make_array(n_x, /fl)
		y_data	= make_array(n_y, /fl)

		readf, lun, x_data
		readf, lun, y_data
		n_xy	= n_x * n_y
		isort	= -1
	end
	; irregular grid, user defined labels
	20: begin
		siscape, lun, x_title, n_xy
		siscape, lun, y_title, n_xy

		x_data	= make_array(n_xy, /fl)
		y_data	= make_array(n_xy, /fl)

		tmp	= make_array(2, /fl)
		for i=0L,n_xy-1 do begin
			readf, lun, tmp
			x_data(i)	= tmp(1)
			y_data(i)	= tmp(0)
		endfor
		isort	= 0
	end
	endcase

	irr	= 1 - abs(isort)
	free_lun, lun
end; read_grid
;-----------------------------------------------------------------------



;- data ----------------------------------------------------------------
pro read_data, data_file, idata, log10=log10, zrange=zrange, scale=scale
	common	data,	z_data, z_title, u_title
	common	grid
	;columns	= 4
	;n_data	= FILE_LINES(data_file)
	result	= QUERY_ASCII(data_file, info)
	if (result ne 1) then begin
		print, 'error in QUERY_ASCII on file:', data_file
		stop
	endif
	n_data	= info.lines
	columns	= info.words / n_xy
	;columns	= info.words / n_data
	print, 'columns in data file: ', columns

;	if (n_data ne n_xy) then begin
;		print, 'resolution of grid file (',n_xy,') does not match data resolution (',n_data,')'
;		stop
;	endif

	f_scale	= 1.0
	if (keyword_set(scale)) then f_scale = scale

	data	= make_array(n_xy, columns, /fl)
	tmp	= make_array(columns, /fl)
	pi2	= 2*!PI

;	z_title_strs	= ['density',$
;			   'el. temperature',$
;			   'mach number',$
;			   'bfield strength']
;	u_title_strs	= ['n_e [10!E19!N m!E3!N]',$
;			   'T_e [eV]',$
;			   '',$
;			   'B [T]']
;	u_title	= u_title_strs(idata)
;	z_title	= z_title_strs(idata)
	u_title	= ''
	z_title	= ''

	openr, lun, data_file, /get_lun
	i=0L
	str=' '
	while not eof(lun) do begin
	;for i=0L,n_data-1 do begin
		;readf, lun, tmp
		;data(i,0) = tmp(0)
		;data(i,0) = tmp(0) * 10
		;data(i,0) = tmp(0) * 22.5
		;data(i,1) = tmp(1)
		;data(i,2) = tmp(2)
		;data(i,2) = tmp(2) * 1000
		;data(i,3) = tmp(3)
		readf, lun, str
		if (strmid(str,0,1) ne '#') then begin
			reads, str, tmp
			data(i,*) = tmp(*) * f_scale
			i	= i + 1
		endif
		if (i eq n_xy) then break
	endwhile
	if (i ne n_xy) then begin
		print, 'resolution of grid file (',n_xy,') does not match data resolution: ', i
		stop
	endif
	;endfor
	free_lun, lun

	if (keyword_set(log10)) then begin
	for i=0L,n_xy-1 do begin
		d	= data(i,idata)
		if (d gt 0.0) then begin
			data(i,idata)	= alog10(d)
		endif else begin
			data(i,idata)	= -100
		endelse

		d	= data(i,idata)
		if (keyword_set(zrange)) then begin
			if (d lt zrange(0)) then data(i,idata) = zrange(0)
			if (d gt zrange(1)) then data(i,idata) = zrange(1)
		endif
	endfor
	endif

	case isort of
	; irregular grid
	0: begin
		z_data	= data(*,idata)
	end
	; regular grid (in columns)
	1: begin
		z_data	= make_array(n_x, n_y, /fl)
		for i=0L,n_x-1 do begin
		for j=0L,n_y-1 do begin
			z_data(i,j)	= data(i*n_y + j,idata)
		endfor
		endfor
	end
	; regular grid (in rows)
	-1: begin
		z_data	= make_array(n_x, n_y, /fl)
		for i=0L,n_x-1 do begin
		for j=0L,n_y-1 do begin
			z_data(i,j)	= data(i + j*n_x,idata)
		endfor
		endfor
	end
	endcase
end; read_data
;-----------------------------------------------------------------------



;- device --------------------------------------------------------------
pro open_device, ps_plot=ps_plot, ct=ct, white0=white0, xsize=xsize, ysize=ysize
	common	cb_plot,	char_thick, char_size, $
				XL_Margin, XR_Margin, XL_Margin_CBar, XR_Margin_CBar

if (keyword_set(ps_plot)) then begin
	set_plot, 'PS'
	ps_file_	= ps_plot

	if (not keyword_set(xsize)) then xsize	= 16
	if (not keyword_set(ysize)) then ysize	= 12
	;device, file=ps_file_, /ENCAPSUL, /COLOR, BITS=8, XSIZE=20, YSIZE=12
	;device, file=ps_file_, /ENCAPSUL, /COLOR, BITS=8, XSIZE=16, YSIZE=12
	device, file=ps_file_, /ENCAPSUL, /COLOR, BITS=8, XSIZE=xsize, YSIZE=ysize
	XL_Margin	= 8
	XR_Margin	= 14
	;XL_Margin_CBar	= 60;		for charsize=1.4, charthick=3.0, xsize=20
	;XL_Margin_CBar	= 48;		for charsize=1.4, charthick=3.0, xsize=16
	XL_Margin_CBar	= xsize*3;	for charsize=1.4, charthick=3.0
	;XL_Margin_CBar	= 85
	XR_Margin_CBar	= 2
	char_size	= 1.4
	char_thick	= 3.0
endif else begin
	set_plot, 'X'
	device, retain=2, decomposed=0
	base	= WIDGET_BASE (/ROW, TITLE='plot window', /BASE_ALIGN_CENTER)
	plot_window	= WIDGET_DRAW (base, XSIZE=800, YSIZE=600)
	WIDGET_CONTROL, base, /REALIZE
	XL_Margin	= 8
	XR_Margin	= 23
	XL_Margin_CBar	= 127
	XR_Margin_CBar	= 3
	char_size	= 1.0
	char_thick	= 1.0
endelse
	;loadct, 39, file='/home/frerichs/work/zTemplates/color_tables/CT_emc3_display'
	;loadct, 45, file='/home/frerichs/work/zTemplates/color_tables/CT_emc3_display'	; depostion table
	ct0=39
	if (keyword_set(ct)) then ct0=ct
	;loadct, ct0, file='/home/frerichs/work/zTemplates/color_tables/CT_emc3_display'
	loadct, ct0
	tvlct, ct_red, ct_green, ct_blue, /GET

	if (keyword_set(white0)) then begin
		ct_red(1)	= 255
		ct_green(1)	= 255
		ct_blue(1)	= 255
	endif
	;ct_red(0)	= 0
	;ct_green(0)	= 0
	;ct_blue(0)	= 0
	;ct_red(255)	= 255
	;ct_green(255)	= 255
	;ct_blue(255)	= 255
	tvlct, ct_red, ct_green, ct_blue
end; open_device
pro close_device, ps_plot=ps_plot
	if (keyword_set(ps_plot)) then device, /close
end; close_device
;-----------------------------------------------------------------------



;- poincare ------------------------------------------------------------
pro plot_poincare, poincare_file
	common	grid
	print, 'overlaying poincare plot'
	columns	= 5
	n_data	= FILE_LINES(poincare_file)

	tmp	= make_array(columns, /fl)
	p_data	= make_array(n_data, 2, /fl)

	openr, lun, poincare_file, /get_lun
	for i=0L,n_data-1 do begin
		readf, lun, tmp
		
		switch grid_id of
		2:
		4: begin
			p_data(i,0)	= 360.0 - tmp(2)
			p_data(i,1)	= tmp(3)
			break
		end
		7:
		8: begin
			p_data(i,0)	= tmp(2)
			p_data(i,1)	= tmp(4)
			break
		end
		else:	p_data(i,*)	= tmp(0:1)
		endswitch
	endfor

	oplot, p_data(*,0), p_data(*,1), psym=3, color=0
	free_lun, lun
end; plot_poincare
;-----------------------------------------------------------------------



;-----------------------------------------------------------------------
; overlay additional data				2010-05-27
;
; generalization of plot_poincare
pro overlay_data, data_files, plot_style=plot_style
	common  grid

	print, 'plotting additional data as overlay ...'

	n	= (size(data_files))[1]
	for i=0,n-1 do begin
		file_str	= data_files(i)

		; scan data file for data lines
		str		= ''
		nd1		= FILE_LINES(file_str)
		nd2		= 0L
		openr, lun, file_str, /get_lun
		for j=1,nd1 do begin
			readf, lun, str
			if (strmid(str,0,1) ne '#') then nd2 = nd2 + 1
		endfor
		free_lun, lun

		; read data
		xy_data		= make_array(2, nd2, /fl)
		nd2		= 0L
		x		= 0.0
		y		= 0.0
		openr, lun, file_str, /get_lun
		for j=1,nd1 do begin
			readf, lun, str
			if (strmid(str,0,1) ne '#') then begin
				reads, str, x, y
				xy_data(0,nd2) = x 
				xy_data(1,nd2) = y 
				nd2 = nd2 + 1
			endif
		endfor
		free_lun, lun


		; select plot style
		if not keyword_set(plot_style) then begin
			style	= 0
		endif else begin
			style	= plot_style(i)
		endelse

		; plot data
		case style of
		0:	oplot, xy_data(0,*), xy_data(1,*), thick=3.0
		1:	polyfill, xym_data(0,*), xy_data(1,*), color=fill_color, noclip=0
		2:	oplot, xy_data(0,*), xy_data(1,*), thick=3.0, color=0
		2:	oplot, xy_data(0,*), xy_data(1,*), thick=3.0, color=0
		3:	oplot, xy_data(0,*), xy_data(1,*), thick=3.0, color=251, psym=1
		4:	oplot, xy_data(0,*), xy_data(1,*), thick=3.0, color=252, psym=1
		5:	oplot, xy_data(0,*), xy_data(1,*), thick=3.0, color=253, psym=1
		6:	oplot, xy_data(0,*), xy_data(1,*), thick=3.0, color=254, psym=1
		7:	oplot, xy_data(0,*), xy_data(1,*), thick=10.0, color=255
		8:	oplot, xy_data(0,*), xy_data(1,*), color=252, psym=3
		9:	oplot, xy_data(0,*), xy_data(1,*), thick=10.0, color=251
		10:	oplot, xy_data(0,*), xy_data(1,*), thick=10.0, color=0, linestyle=2
		100:	oplot, xy_data(0,*), xy_data(1,*), thick=3.0, color=0, psym=6, symsize=2
		endcase
	endfor
end; overlay_data
;-----------------------------------------------------------------------



;- boundary ------------------------------------------------------------
pro plot_boundary, boundary_file
	n	= (size(boundary_file))[1]

	for i=0,n-1 do begin
		file_str	= boundary_file(i)
		n_data		= FILE_LINES(file_str)
		lim_data	= make_array(2, n_data, /fl)
		openr, lun, file_str, /get_lun
		readf, lun, lim_data
		free_lun, lun
		;lim_data(0,*)	= 360.0 - lim_data(0,*)
		oplot, lim_data(0,*), lim_data(1,*), psym=0, color=255
		;polyfill, lim_data(0,*), lim_data(1,*), color=255, noclip=0
	endfor
end; plot_boundary
;-----------------------------------------------------------------------



;- overly contour of connection length ---------------------------------
pro plot_lc_contour, lc_file, lc_data
	common	grid
	common	cb_plot

	read_grid, 'grid_hr_yx.dat'
	columns	= 6
	n_data	= FILE_LINES(lc_file)

	if (n_data ne n_xy) then begin
		print, 'resolution of grid file (',n_xy,') does not match data resolution (',n_data,')'
		stop
	endif

	data	= make_array(n_data, /fl)
	tmp	= make_array(columns, /fl)
	pi2	= 2*!PI

	openr, lun, lc_file, /get_lun
	for i=0L,n_data-1 do begin
		readf, lun, tmp
		lc1	= abs(tmp(0))
		lc2	= abs(tmp(1))
		data(i) = (lc1 + lc2) /1.e2	; total connection length [m]

		np1	= abs(tmp(2))
		np2	= abs(tmp(3))
		data(i) = (np1 + np2)/pi2		; total poloidal turns

	endfor
	free_lun, lun

; sort data according to grid
	case isort of
	; irregular grid
	0: begin
		lc_data	= data
	end
	; regular grid (in columns)
	1: begin
		lc_data	= make_array(n_x, n_y, /fl)
		for i=0L,n_x-1 do begin
		for j=0L,n_y-1 do begin
			lc_data(i,j)	= data(i*n_y + j)
		endfor
		endfor
	end
	; regular grid (in rows)
	-1: begin
		lc_data	= make_array(n_x, n_y, /fl)
		for i=0L,n_x-1 do begin
		for j=0L,n_y-1 do begin
			lc_data(i,j)	= data(i + j*n_x)
		endfor
		endfor
	end
	endcase

	;clevels=[50,100,200]
	clevels=[1,2,3]
	contour, lc_data, x_data, y_data, irregular=irr, xstyle=5, ystyle=5, $
		 xmargin=[XL_Margin,XR_Margin], levels=clevels, /noerase, color=255, c_thick=1.0, $
		 charsize=char_size, charthick=char_thick
end; plot_lc_contour
;-----------------------------------------------------------------------



;-----------------------------------------------------------------------
pro write_info
end; write_info
;-----------------------------------------------------------------------


;- plot_data -----------------------------------------------------------
; visualization tool for 2D data
;--------------------------
;- required parameters
;
; grid_file:	file with definitions of grid nodes
; data_file:	file with data values on grid nodes (one row per node)
; idata:	select column # for plotting
;
;--------------------------
;- optional parameters
;
; white0:	plot lowest value in white (useful for printable output)
;
; clevels:	array of contour levels to be plotted
;
; xselect,
; yselect:	extract 1D profiles at given x (y) position (only for regular grids, i.e. grid types 3,4,41,5,50,6,8)
; xfile, yfile:	optional output file name for 1D profiles
;
;- plot_data -----------------------------------------------------------
pro plot_data, grid_file, data_file, idata, zrange=zrange, ps_plot=ps_plot, poincare=poincare_file, boundary=boundary_file, label=label, lpos=lpos, xrange=xrange, yrange=yrange, clevels=clevels, c_labels=c_labels, nlevels=nlevels, title=title, $
		plot_add_data=plot_add_data, add_data_style=add_data_style, log10=log10, scale=scale, $
		xtitle=xtitle, ytitle=ytitle, utitle=utitle, ct=ct, white0=white0, $
		xselect=xselect, yselect=yselect, xfile=xfile, yfile=yfile, xsize=xsize, ysize=ysize
	common	grid
	common	data
	common	cb_plot

	if (idata lt 0)	then begin
		write_info
		return
	endif

	read_grid, grid_file

	read_data, data_file, idata, log10=log10, zrange=zrange, scale=scale

	open_device, ps_plot=ps_plot, ct=ct, white0=white0, xsize=xsize, ysize=ysize

	if (keyword_set(zrange)) then begin
		z_min	= zrange(0)
		z_max	= zrange(1)
	endif else begin
		z_min	= min(z_data)
		z_max	= max(z_data)
	endelse

	; number of color levels
	if (keyword_set(nlevels)) then begin
		n_clevels	= nlevels + 1
	endif else begin
		n_clevels	=  64
	endelse

	n_colors	= 253
	color_min	= 1
	levels	= z_min + 1.0*findgen(n_clevels)/(n_clevels-1)*(z_max-z_min)

	colors	= 1.0*findgen(n_clevels)/(n_clevels-1) * n_colors + color_min

	z_title	= ''
	u_title = ''
	;z_title	= 'Magnetic topology (Textor shot 102168)'
	;z_title	= 'Magnetic topology (DIII-D shot 125992)'
	;z_title	= '   Magnetic topology (DIII-D n=2 scenario)'
	;z_title	= 'Magnetic footprint (DIII-D shot 125992)'
	;z_title	= 'backward direction'
	;z_title	= 'forward direction'
	;z_title	= 'L_c at toroidal plane'
	;z_title	= 'Magnetic footprint (DIII-D shot 122342)'
	;z_title	= 'Magnetic footprint on ALT-II limiter (shot 102087)'
	;z_title	= 'Magnetic footprint on DED limiter (shot 102087)'
	;z_title	= 'Magnetic topology (shot 102087)'
	;z_title	= 'Magnetic topology (shot 96896 + 6/2)'
	;y_title	= 'Z on wall at HFS [cm]'
	;y_title	= 'position on wall [cm]'
	;y_title	= 'poloidal angle [deg]'
	;z_title	= 'Magnetic topology (shot 105377, I_DED = 2.3 kA)'
	;z_title	= 'Magnetic topology (shot 105377, I_DED = 2.3 kA)'
	;z_title	= 'H!D!4a!X!N Emission'
	;u_title = 'P_rad [a.u.]'
	;u_title = 'n!De!N [10!U19!N m!U-3!N]'
	;u_title = '        [kW m!U-3!N m]'
	;u_title = 'P!Drad!N [kW m!U-2!N]'
	;u_title = 'P!Drad!N [kW m!U-3!N]'
	;z_title		= 'Ionization source'
	;u_title		= 'S!Dp!N [kA m!U-3!N]'

	if (keyword_set(title)) then z_title = title
	if not keyword_set(xtitle) then xtitle = x_title
	if not keyword_set(ytitle) then ytitle = y_title
	if not keyword_set(utitle) then utitle = u_title

	; plot actual data
	contour, z_data, x_data, y_data, irregular=irr, /fill, xstyle=1, ystyle=1, $
                 xmargin=[XL_Margin,XR_Margin], levels=levels, c_colors=colors, $
		 xtitle=xtitle, ytitle=ytitle, title=z_title, $
		 xrange=xrange, yrange=yrange, $
		 charsize=char_size, charthick=char_thick ;, /isotropic ??????

	; contour lines
	if (keyword_set(clevels))	then begin
		contour, z_data, x_data, y_data, irregular=irr, xstyle=5, ystyle=5, $
			 xmargin=[XL_Margin,XR_Margin], levels=clevels, /noerase, color=0, c_thick=2.0, $
			 charsize=char_size, charthick=char_thick, c_labels=c_labels
;		contour, z_data, x_data, y_data, irregular=irr, xstyle=5, ystyle=5, $
;			 xmargin=[XL_Margin,XR_Margin], levels=clevels, /noerase, color=255, c_thick=2.0, $
;			 charsize=char_size, charthick=char_thick, c_labels=c_labels, c_linestyle=2
	endif

	;plot_lc_contour, 'LAMINAR_PLOT'

	;oplot, [120,240], [47.5, 47.5], color=255, thick=2, linestyle=2

	; ALTII-corners
	;oplot, [-22.5,157.5], [332.24, 332.24], color=255, thick=8, linestyle=2
	;oplot, [-22.5,157.5], [297.76, 297.76], color=255, thick=8, linestyle=2

	if (keyword_set(poincare_file)) then plot_poincare, poincare_file
	if (keyword_set(boundary_file)) then plot_boundary,  boundary_file
	if (keyword_set(label))		then xyouts, lpos(0), lpos(1), label, charsize=1.5, charthick=2.0
	if (keyword_set(plot_add_data))	then overlay_data,  plot_add_data, plot_style=add_data_style

	;xyouts, 315, 1.04, 'ALT-II target', charsize=1.5, charthick=2.0, orientation=90
	;xyouts, 320, 1.02, 'ALT-II', charsize=1.2, charthick=2.0, orientation=90

	; plot color bar
	bar_dummy	= transpose(levels#[1,1])
	colBarStr	= z_title

	contour, bar_dummy, [0,1], levels, /fill, levels=levels, $
		 ystyle=1, xstyle=4, xmargin=[XL_Margin_CBar,XR_Margin_Cbar], /noerase, $
		 ytitle=utitle, c_colors=colors, $
		 charsize=char_size, charthick=char_thick


	close_device, ps_plot=ps_plot


	if (isort ne 0) then begin
	; extract 1D profile at given x
	if (keyword_set(xselect)) then begin
		if (xselect ge 0  and  xselect lt n_x) then begin
			output_file	= 'xprofile.txt'
			if (keyword_set(xfile)) then output_file	= xfile
			openw, lun, output_file, /get_lun
			printf, lun, '# 1D profile at x = ', x_data(xselect)
			for iy=0,n_y-1 do begin
				printf, lun, y_data(iy), z_data(xselect, iy)
			endfor
			free_lun, lun
		endif
	endif

	; extract 1D profile at given y
	if (keyword_set(yselect)) then begin
		if (yselect ge 0  and  yselect lt n_y) then begin
			output_file	= 'yprofile.txt'
			if (keyword_set(yfile)) then output_file	= yfile
			openw, lun, output_file, /get_lun
			printf, lun, '# 1D profile at y = ', y_data(yselect)
			for ix=0,n_x-1 do begin
				printf, lun, x_data(ix), z_data(ix, yselect)
			endfor
			free_lun, lun
		endif
	endif
	endif
end; plot_data
;-----------------------------------------------------------------------
