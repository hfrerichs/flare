;- grid ----------------------------------------------------------------
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
		x_title	= ''
		y_title	= ''
		tmp	= make_array(3, /fl)
		for i=0L,n_xyz-1 do begin
			readf, lun, tmp
			x_data(i)	= tmp(0)
			y_data(i)	= tmp(1)
			;y_data(i)	= tmp(2)
		endfor
		n_xy	= n_xyz
		isort	= 0
	end
	1: begin
		iscape, lun, n_RZ
		rscape, lun, phi

		x_data	= make_array(n_RZ, /fl)
		y_data	= make_array(n_RZ, /fl)
		x_title	= 'Major Radius [cm]'
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
		x_title	= 'Poloidal Angle [deg]'
		y_title	= 'Minor Radius [cm]'

		tmp	= make_array(2, /fl)
		for i=0L,n_rt-1 do begin
			readf, lun, tmp
			;x_data(i)	= tmp(1)
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
		x_title	= 'Poloidal Angle [deg]'
		y_title	= 'Minor Radius [cm]'

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
	41: begin;	same as 4 but with transposed x,y coordinates for plotting
		iscape, lun, n_t
		iscape, lun, n_r

		y_data	= make_array(n_t, /fl)
		x_data	= make_array(n_r, /fl)
		y_title	= 'Poloidal Angle [deg]'
		x_title	= 'Minor Radius [cm]'

		rscape, lun, r_center
		rscape, lun, phi

		readf, lun, x_data
		readf, lun, y_data
		y_data	= 360.0 - y_data
		n_y	= n_t
		n_x	= n_r
		n_xy	= n_x * n_y
		isort	= -1
	end
	5: begin
		iscape, lun, n_p
		iscape, lun, n_t

		x_data	= make_array(n_p, /fl)
		y_data	= make_array(n_t, /fl)
		x_title	= 'Toroidal Angle [deg]'
		y_title	= 'Poloidal Angle [deg]'

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
		x_title	= 'Toroidal Angle [deg]'
		y_title	= 'Length [cm]'

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

		x_data	= make_array(n, /fl)
		y_data	= make_array(n, /fl)
		x_title	= 'Poloidal Angle [deg]'
		y_title	= '!4w!3!IN!N (i.e. radial coordinate)'
		;y_title	= 'normalized flux !4w!3!IN!N (radial coord.)'

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
		x_title	= 'Poloidal Angle [deg]'
		y_title	= '!4w!3!IN!N (i.e. radial coordinate)'
		;y_title	= 'normalized flux !4w!3!IN!N (radial coord.)'

		readf, lun, x_data
		readf, lun, y_data
		;x_data	= 360.0 - x_data
		n_x	= n_theta
		n_y	= n_psi
		n_xy	= n_x * n_y
		isort	= -1
	end
	endcase

        if (keyword_set(xoffset) and (isort ne 0)) then begin
		x_data(*)	= x_data(*) + xoffset
	endif
        if (keyword_set(yoffset) and (isort ne 0)) then begin
		y_data(*)	= y_data(*) + yoffset
	endif

	irr	= 1 - abs(isort)
	x_min	= min(x_data)
	x_max	= max(x_data)
	y_min	= min(y_data)
	y_max	= max(y_data)
	free_lun, lun
end; read_grid
;-----------------------------------------------------------------------



;- data ----------------------------------------------------------------
pro read_data, data_file, idata, log10=log10, zrange=zrange
	common	data,	z_data, z_title, u_title
	common	grid


	z_title_strs	= ['Connection Length',$
			   'shortest connection length',$
			   'effective connection length',$
			   'Poloidal Turns',$
			   'shortest poloidal turns',$
			   '- minimum psi',$
			   '- average psi',$
			   'Lc in backward direction',$
			   'Lc in forward direction', $
			   'Lc in shortest direction', $
			   'Penetration depth',$
			   'Wall index',$
			   'Wall index',$
			   'Wall index',$
			   'Pol. locat. of minPsiN',$
			   'Tor. locat. of minPsiN',$
			   'Radial distance to LCFS']
	u_title_strs	= ['Connection Length [m]',$
			   'Shortest Length [m]',$
			   'Length [m]',$
			   'Poloidal Turns',$
			   'Poloidal Turns',$
			   'psiN',$
			   'psiN',$
			   'Lc- [m]',$
			   'Lc+ [m]',$
			   'Lc [p.t.]',$
			   'Penetration depth',$
			   '',$
			   '',$
			   '',$
			   'Pol. locat. of minPsiN',$
			   'Tor. locat. of minPsiN',$
			   'Radial distance to LCFS']


	; check data file for the correct # of lines ..................................
	n_data	= FILE_LINES(data_file)
	if (n_data ne n_xy) then begin
		print, 'resolution of grid file (',n_xy,') does not match data resolution (',n_data,')'
		stop
	endif
	;..............................................................................


	; get configuration of additional data ........................................
	info_file	= data_file+'.info'
	nout		= 1
	iout		= make_array(3, /int)
	iout(*)		= 0
	iout(0)		= 1
	openr, lun, info_file, /get_lun, error=err
	if (err eq 0) then begin
	readf, lun, nout
	j = 0
	a = ''
	for i = 1, nout do begin
		readf, lun, j, a, format='(i4,a32)'
		iout(i-1) = j
	endfor
	free_lun, lun
	endif
	;columns	= 6
	columns	= 5 + nout
	n_vars	= 10 + 4 + 2 + 1
	;..............................................................................
	

	; setup data arrays + additional stuff ........................................
	data	= make_array(n_data, n_vars, /fl)
	tmp	= make_array(columns, /fl)
	pi2	= 2*!PI

	u_title	= u_title_strs(idata)
	z_title	= z_title_strs(idata)
	;..............................................................................
	

	; read data ...................................................................
	openr, lun, data_file, /get_lun
	for i=0L,n_data-1 do begin
		;print, i
		readf, lun, tmp
		;data(i,*) = tmp
		lc1	= abs(tmp(0))
		lc2	= abs(tmp(1))
		data(i,0) = (lc1 + lc2) /1.e2		; total connection length [m]
		;data(i,0) = alog10(data(i,0))
		;data(i,0) = lc2/1.e2
		data(i,1) = (min([lc1, lc2])) /1.e2	; shortest connection length [m]
		;data(i,2) = lc1*lc2 / (lc1 + lc2) /1.e2	; effective connection length [m]
		;data(i,2) = 2 * lc1 / (lc1 + lc2) - 1	; effective connection length
		data(i,2) = (lc2-lc1) / (lc1 + lc2)	; effective connection length
;		if (lc1 gt lc2) then begin
;			data(i,2)	= lc1 /1.e2
;		end else begin
;			data(i,2)	= -lc2 /1.e2
;		endelse


		np1	= abs(tmp(2))
		np2	= abs(tmp(3))
		data(i,3) = (np1 + np2)			; total poloidal turns
		data(i,4) = min([np1, np2])		; shortest poloidal turns
		;data(i,3) = (np1 + np2)/pi2		; total poloidal turns
		;data(i,4) = min([np1, np2])/pi2		; shortest poloidal turns
		data(i,5) = tmp(4)			; minimum psiN
		data(i,10)	= 1.0 - data(i,5)	; 1 - min(psiN)
;		if (data(i,10) le 0.0) then begin
;			data(i,10)	= -99
;		end else begin
;			data(i,10)	= alog10(data(i,10))
;		endelse

		data(i,6) = tmp(5)			; average psiN
		data(i,6) = 0.0			; average psiN

		data(i,7)	= np1
		data(i,8)	= np2
		;data(i,7)	= lc1/1.e2		; Lc in backward direction [m]
		;data(i,8)	= lc2/1.e2		; Lc in forward direction [m]

		if (np1 gt np2) then begin
			data(i,9)	= - data(i,4)
		end else begin
			data(i,9)	=   data(i,4)
		endelse

		data(i,11)	= 0.0
		data(i,12)	= 0.0
		;data(i,11)	= tmp(6)
		;data(i,12)	= tmp(7)

		for j=0,nout-1 do begin
		case iout(j) of
		0: begin
		end
		1: begin
			data(i,6) = tmp(5+j)			; average psiN
		end
		2: begin
			data(i,11)= tmp(5+j)			; Limiting surface # (neg. direction)
		end
		4: begin
			data(i,12)= tmp(5+j)			; Limiting surface # (pos. direction)
		end
		8: begin
			data(i,14)= tmp(5+j)			; Toroidal location of deepest penetration
		end
		16: begin
			data(i,15)= tmp(5+j)			; Poloidal location of deepest penetration
		end
		512: begin
			data(i,16)= tmp(5+j)			; radial distance to LCFS
		end
		endcase
		endfor

		data(i,13)	= data(i,12)
		;if (data(i,13) eq -1) then data(i,13) = data(i,11)
		;if (data(i,8) le 1.0) then data(i,13) = data(i,11)
		if (data(i,8) le data(i,7)) then data(i,13) = data(i,11)
	endfor
	free_lun, lun
	;..............................................................................


	; plot logarithmic data if requested ..........................................
	if (keyword_set(log10)) then begin
	for i=0L,n_data-1 do begin
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
	;..............................................................................


	; select requested data (and build data structure) ............................
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
	;..............................................................................
end; read_data
;-----------------------------------------------------------------------



;- device --------------------------------------------------------------
pro set_ct, ct_red, ct_green, ct_blue
	A		= make_array(15, /fl)
	B		= make_array(15, /fl)
	r_red		= make_array(65, /fl)
	r_green		= make_array(65, /fl)
	r_blue		= make_array(65, /fl)
	for i=1,15 do begin
		A(i-1)	= i * 1.0 / 16.0
		B(i-1)	= (16-i) * 1.0 / 16.0
	endfor
	r_red(1:24)	= 0.0
	r_red(25:39)	= A
	r_red(40:56)	= 1.0
	r_red(57:63)	= B(0:6)
	r_red(64)	= 0.5
	
	r_green(1:8)	= 0.0
	r_green(9:23)	= A
	r_green(24:40)	= 1.0
	r_green(41:55)	= B
	r_green(56:63)	= 0.0
	r_green(64)	= 0.0

	r_blue(1:7)	= A(8:14)
	r_blue(8:24)	= 1.0
	r_blue(25:39)	= B
	r_blue(40:63)	= 0.0
	r_blue(64)	= 0.0

	ct_red(1:64)	= fix(255*r_red(1:64))
	ct_green(1:64)	= fix(255*r_green(1:64))
	ct_blue(1:64)	= fix(255*r_blue(1:64))
end; set_ct

pro open_device, ps_plot=ps_plot, ps_geo=ps_geo, cb_center=cb_center, cb_cutoff=cb_cutoff, color_table=color_table, reverse_color=reverse_color
	common	cb_plot,	char_thick, char_size, $
				XL_Margin, XR_Margin, $
				XL_Margin_CBar, XR_Margin_CBar, XU_Margin_CBar, XD_Margin_CBar

	XU_Margin_CBar	= 2;	default values
	XD_Margin_CBar	= 4
if (keyword_set(ps_plot)) then begin
	set_plot, 'PS'
	ps_file_	= ps_plot

	xsize	= 20
	ysize	= 12
	if (keyword_set(ps_geo)) then begin
		xsize	= ps_geo(0)
		ysize	= ps_geo(1)
	endif
	device, file=ps_file_, /encapsul, /color, bits=8, xsize=xsize, ysize=ysize

	;device, file=ps_file_, /ENCAPSUL, /COLOR, BITS=8, XSIZE=16, YSIZE=24
	;device, file=ps_file_, /ENCAPSUL, /COLOR, BITS=8, XSIZE=20, YSIZE=24
	;device, file=ps_file_, /ENCAPSUL, /COLOR, BITS=8, XSIZE=20, YSIZE=12
	;device, file=ps_file_, /ENCAPSUL, /COLOR, BITS=8, XSIZE=16, YSIZE=12

	if (keyword_set(cb_center)) then begin
		XL_Margin	= 8
		XR_Margin	= 3
		;XL_Margin_CBar	= 30;		16x24 with frame
		;XR_Margin_CBar	= 20
		;XU_Margin_CBar	= 10
		;XD_Margin_CBar	= 18
		;XL_Margin_CBar	= 26;		16x24 without frame
		;XR_Margin_CBar	= 24
		;XU_Margin_CBar	= 13.5
		;XD_Margin_CBar	= 15.5
		XL_Margin_CBar	= 26;		16x28 without frame
		XR_Margin_CBar	= 24
		XU_Margin_CBar	= 14.5
		XD_Margin_CBar	= 16.5
	endif else begin
		XL_Margin	= 8
		XR_Margin	= 14
		;XL_Margin_CBar	= 60;		20 x 12
		;XR_Margin_CBar	= 2
		;XL_Margin_CBar	= 48;		16 x 12
		;XR_Margin_CBar	= 2
		XL_Margin_CBar	= 3*xsize
		XR_Margin_CBar	= 2

		;XL_Margin_CBar	= 85
	endelse
	char_size	= 1.4
	char_thick	= 3.0
endif else begin
	set_plot, 'X'
	device, retain=2, decomposed=0
	base	= WIDGET_BASE (/ROW, TITLE='plot window', /BASE_ALIGN_CENTER)
	plot_window	= WIDGET_DRAW (base, XSIZE=800, YSIZE=600)
	WIDGET_CONTROL, base, /REALIZE

	if (keyword_set(cb_center)) then begin
		XL_Margin	= 8
		XR_Margin	= 8
		XL_Margin_CBar	= 55
		XR_Margin_CBar	= 65
	endif else begin
		XL_Margin	= 8
		XR_Margin	= 23
		XL_Margin_CBar	= 127
		XR_Margin_CBar	= 3
	endelse
	char_size	= 1.0
	char_thick	= 1.0
endelse
	if (keyword_set(color_table)) then begin
		;loadct, 45, file='/home/frerichs/work/zTemplates/color_tables/CT_emc3_display'
		;loadct, color_table, file='/home/frerichs/work/zTemplates/color_tables/CT_emc3_display'
		;loadct, color_table, file='/home/frerichs/sandbox/MagFTC/templates/plot/COLOR_TABLES'
		loadct, color_table
		tvlct, ct_red, ct_green, ct_blue, /GET

		; compress color table
		n		= 64
		rgb_data	= make_array(3, n+1, /fl)
		for i=0,n do begin
			rgb_data(0,i) = ct_red(fix(round(i*255.0/n)))
			rgb_data(1,i) = ct_green(fix(round(i*255.0/n)))
			rgb_data(2,i) = ct_blue(fix(round(i*255.0/n)))
		endfor
		for i=0,n do begin
			ct_red(i)	= rgb_data(0,i)
			ct_green(i)	= rgb_data(1,i)
			ct_blue(i)	= rgb_data(2,i)
		endfor
	endif else begin
		tvlct, ct_red, ct_green, ct_blue, /GET
		set_ct, ct_red, ct_green, ct_blue
	endelse

	ct_red(0)	= 0
	ct_green(0)	= 0
	ct_blue(0)	= 0

	if (keyword_set(cb_cutoff)) then begin
	; lower cut-off
	if (cb_cutoff(0) ge 0 and cb_cutoff(0) le 255) then begin
		;ct_red(1)	= 255	; lower limit = white for printing
		;ct_green(1)	= 255
		;ct_blue(1)	= 255
		ct_red  ( 64)	= ct_red  (cb_cutoff(0))
		ct_green( 64)	= ct_green(cb_cutoff(0))
		ct_blue ( 64)	= ct_blue (cb_cutoff(0))
	endif

	; uppter cut-off
	if (cb_cutoff(1) ge 0 and cb_cutoff(1) le 255) then begin
		ct_red  ( 64)	= ct_red  (cb_cutoff(1))
		ct_green( 64)	= ct_green(cb_cutoff(1))
		ct_blue ( 64)	= ct_blue (cb_cutoff(1))
	endif
	endif

	; reverse
	if (keyword_set(reverse_color)) then begin
	for i=1,32 do begin
		ir	= ct_red  (i)
		ig	= ct_green(i)
		ib	= ct_blue (i)
		ct_red  (i)	= ct_red  (65-i)
		ct_green(i)	= ct_green(65-i)
		ct_blue (i)	= ct_blue (65-i)
		ct_red  (65-i)	= ir
		ct_green(65-i)	= ig
		ct_blue (65-i)	= ib
	endfor
	endif


	ct_red  (250)	= 220
	ct_green(250)	= 220
	ct_blue (250)	= 220
; temporary for TEXTOR fluxtubes
	ct_red  (251)	= 255
	ct_green(251)	=   0
	ct_blue (251)	=   0
	ct_red  (252)	=   0
	ct_green(252)	= 255
	ct_blue (252)	=   0
	ct_red  (253)	=   0
	ct_green(253)	= 255
	ct_blue (253)	= 255
	;ct_red  (254)	= 255
	;ct_green(254)	=   0
	;ct_blue (254)	= 255

	ct_red(254)	= 191
	ct_green(254)	= 191
	ct_blue(254)	= 191
	ct_red(255)	= 255
	ct_green(255)	= 255
	ct_blue(255)	= 255
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



;- limiter -------------------------------------------------------------
pro plot_limiter, limiter_file, fill_color=fill_color, lim_style=lim_style
	n	= (size(limiter_file))[1]

	if (not keyword_set(fill_color)) then fill_color=255
	for i=0,n-1 do begin
		str		= ' '
		file_str	= limiter_file(i)
		; scan file for data lines
		n_data		= 0L
		openr, lun, file_str, /get_lun
		while not eof(lun) do begin
			readf, lun, str
			if (strmid(str,0,1) ne '#') then n_data = n_data + 1
		endwhile
		free_lun, lun

		; setup array for xy data
		;n_data		= FILE_LINES(file_str)-1
		lim_data	= make_array(2, n_data, /fl)
		dline		= make_array(2, /fl)

		openr, lun, file_str, /get_lun
		i_data		= 0L
		while not eof(lun) do begin
			readf, lun, str
			if (strmid(str,0,1) ne '#') then begin
				reads, str, dline
				lim_data(*,i_data)	= dline
				i_data 			= i_data + 1
			endif
		endwhile
		;readf, lun, lim_data
		free_lun, lun
		print, 'plotting limiter from file ', file_str, 'with ', n_data, ' data lines'

		;lim_data(0,*)	= 360.0 - lim_data(0,*)
		;for j=0,n_data-1 do begin
			;if ((lim_data(0,j) gt 180) and (i ge 2)) then lim_data(0,j) = lim_data(0,j) - 360
		;endfor

		;polyfill, lim_data(0,*), lim_data(1,*), color=fill_color, noclip=0
		;oplot, lim_data(0,*), lim_data(1,*), psym=0, color=0, thick=3.2

		if not keyword_set(lim_style) then begin
			style	= 0
		endif else begin
			style	= lim_style(i)
		endelse

;		temporary adjustment for the generation of full scale connection length plots
;		for DIII-D discharge 132731
                a = 8.0
                if (style eq 101) then a = 4.0
		mX = [-1,  1,  1, -1, -1]
		mX = mX / a
		mY = [-1, -1,  1,  1, -1]
		mY = mY / a
		usersym, mX, mY, /fill

		;lim_data(0,*)	= 360.0 - lim_data(0,*)
		case style of
		0:	oplot, lim_data(0,*), lim_data(1,*), thick=3.0
		1:	polyfill, lim_data(0,*), lim_data(1,*), color=fill_color, noclip=0
		2:	oplot, lim_data(0,*), lim_data(1,*), thick=3.0, color=0
		3:	oplot, lim_data(0,*), lim_data(1,*), thick=3.0, color=251, psym=1
		4:	oplot, lim_data(0,*), lim_data(1,*), thick=3.0, color=252, psym=1
		5:	oplot, lim_data(0,*), lim_data(1,*), thick=3.0, color=253, psym=1
		6:	oplot, lim_data(0,*), lim_data(1,*), thick=3.0, color=254, psym=1
		7:	oplot, lim_data(0,*), lim_data(1,*), thick=10.0, color=255
		8:	oplot, lim_data(0,*), lim_data(1,*), color=252, psym=3
		9:	oplot, lim_data(0,*), lim_data(1,*), thick=10.0, color=251
		10:	oplot, lim_data(0,*), lim_data(1,*), thick=10.0, color=0, linestyle=2

		; point styles
		100:	oplot, lim_data(0,*), lim_data(1,*), color=0, psym=8, thick=1.0
		101:	oplot, lim_data(0,*), lim_data(1,*), color=250, psym=8, thick=1.0 
		endcase

	endfor

	;px_arr	= [120, 140, 140, 120, 120]
	;py_arr	= [-120, -120, -100, -100, -120]
	;px_arr	= [100, 160, 160, 100, 100]
	;py_arr	= [-140, -140, -80, -80, -140]
	;oplot, px_arr, py_arr, color=0, thick=8
end; plot_limiter
;-----------------------------------------------------------------------



;-----------------------------------------------------------------------
pro plot_grid, grid_file
	n_data	= FILE_LINES(grid_file)
	xy  	= make_array(2, n_data, /fl)
	str	= ''
	n_block	= -1L
	x	= 0.0
	y	= 0.0

	openr, lun, grid_file, /get_lun
	while ~EOF(lun) do begin
		readf, lun, str
		if ((str ne '') and (str ne ' ')) then begin
			n_block = n_block + 1
			reads, str, x, y
			xy(0, n_block) = x
			xy(1, n_block) = y
		endif else begin
			oplot, xy(0,0:n_block), xy(1,0:n_block)
			n_block	= -1L
		endelse
	endwhile
	free_lun, lun
end; plot_grid
;-----------------------------------------------------------------------



;-----------------------------------------------------------------------
pro print_help
	print, '========================================================================'
	print, ' plot_lc: IDL tool for visualization of the magnetic field structure'
	print, '========================================================================'
	print, ' how to run plot_lc:'
	print
	print, '	plot_lc, grid_file, data_file'
	print
	print, ' optional arguments:'
	print
	print, "	z		= 'Lc':  plot connection length"
	print, "			= 'Lcs': plot shortest connection length"
	print, ' ...'
	print
end; print_help
;-----------------------------------------------------------------------


;- plot_lc -------------------------------------------------------------
pro plot_lc, grid_file, data_file, idata, $
	help=help, $		; print help if keyword set
	Zdata=Zdata, $		; select data for plotting
	ps_plot=ps_plot, ps_geo=ps_geo, poincare=poincare_file, $
	limiter=limiter_file, lim_style=lim_style, fill_color=fill_color, label=label, lpos=lpos, $
	xrange=xrange, xtitle=xtitle, xticks=xticks, xtickv=xtickv, xtickinterval=xtickinterval, xstyle=xstyle, $
	yrange=yrange, ytitle=ytitle, yticks=yticks, ytickv=ytickv, ytickinterval=ytickinterval, ystyle=ystyle, $
	zrange=zrange, ztitle=ztitle, zticks=zticks,                ztickinterval=ztickinterval, $
	xoffset=xoffset, yoffset=yoffset, $
	clevels=clevels, nlevels=nlevels, title=title,$
	plot_grid=grid_fileX, $
	mark_x=mark_x, mark_y=mark_y, $
	cb_center=cb_center, cb_cutoff=cb_cutoff, no_cb=no_cb, color_table=color_table, reverse_color=reverse_color, $
	log10=log10
	common	grid
	common	data
	common	cb_plot

; check user input
	; print help
	if (keyword_set(help)) then begin
		print_help
		return
	endif
	; check mandatory arguments grid_file and data_file
	if (not keyword_set(grid_file)) then begin
		print, 'error: grid file not defined!'
		print, 'try "plot_lc, /help'
		return
	endif
	if (not keyword_set(data_file)) then begin
		print, 'error: data file not defined!'
		print, 'try "plot_lc, /help'
		return
	endif

	idata = 0
	if (keyword_set(Zdata)) then case Zdata of
		'Lc':  idata=0
		'Lcs': idata=1
		'Lceff': idata=2
		'Lpt': idata=3
		'Lpts': idata=4
		'minPsiN': idata=5
		'1minPsiN': idata=10
		'avPsiN': idata=6
		'Lc-': idata=7
		'Lc+': idata=8
		'W-': idata=11
		'W+': idata=12
		'Wlong': idata=13
		'PolLocminPsiN': idata=14
		'TorLocminPsiN': idata=15
		'RD2LCFS': idata=16
		else: begin
			print, 'error: invalid expression for "Zdata"'
			return
		end
	endcase
	if (idata lt 0)	then begin
		write_info
		return
	endif

	read_grid, grid_file, xoffset=xoffset, yoffset=yoffset

	read_data, data_file, idata, log10=log10, zrange=zrange

	open_device, ps_plot=ps_plot, ps_geo=ps_geo, cb_center=cb_center, cb_cutoff=cb_cutoff, color_table=color_table, reverse_color=reverse_color

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

	;n_colors	= 253
	;n_clevels	= 64
	n_colors	= 64
	color_min	= 1
	levels	= z_min + 1.0*findgen(n_clevels)/(n_clevels-1)*(z_max-z_min)
	print, n_clevels, ' contour levels from ',z_min,' to ',z_max

	;colors	= 1.0*findgen(n_clevels-1)/(n_clevels-2) * (n_colors-1) + color_min
	colors	= make_array(n_clevels, /fl)
	colors(0:n_clevels-2)	= 1.0*findgen(n_clevels-1)/(n_clevels-2) * (n_colors-1) + color_min
	colors(n_clevels-1)	= 255

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
	;z_title	= 'Magnetic footprint on ALT-II limiter (scenario B)'
	;z_title	= 'Magnetic footprint on DED limiter (scenario A)'
	;z_title	= 'Magnetic topology (shot 102087)'
	;z_title	= 'Magnetic topology (shot 96896 + 6/2)'
	;y_title	= 'Z on limiting wall [cm]'
	;y_title	= 'position on wall [cm]'
	;y_title	= 'poloidal angle [deg]'
	;z_title	= 'Magnetic topology (shot 105377, I_DED = 2.3 kA)'
	z_title	= ''
	if (keyword_set(title)) then z_title = title
	if not keyword_set(xtitle) then xtitle = x_title
	if not keyword_set(ytitle) then ytitle = y_title
	if not keyword_set(xstyle) then xstyle = 1
	if not keyword_set(ystyle) then ystyle = 1

	; plot actual data
	contour, z_data, x_data, y_data, irregular=irr, /fill, xstyle=xstyle, ystyle=ystyle, $
                 xmargin=[XL_Margin,XR_Margin], levels=levels, c_colors=colors, $
		 xtitle=xtitle, ytitle=ytitle, title=z_title, $
		 xrange=xrange, xticks=xticks, xtickv=xtickv, xtickinterval=xtickinterval, $
		 yrange=yrange, yticks=yticks, ytickv=ytickv, ytickinterval=ytickinterval, $
		 charsize=char_size, charthick=char_thick

	; contour lines
	if (keyword_set(clevels))	then begin
		contour, z_data, x_data, y_data, irregular=irr, xstyle=5, ystyle=5, $
			 xmargin=[XL_Margin,XR_Margin], levels=clevels, /noerase, color=255, c_thick=4.0, $
		 	xrange=xrange, xticks=xticks, xtickv=xtickv, xtickinterval=xtickinterval, $
		 	yrange=yrange, yticks=yticks, ytickv=ytickv, ytickinterval=ytickinterval, $
			 charsize=char_size, charthick=char_thick
	endif

	; mark y-positions
	if (keyword_set(mark_y))	then begin
	;oplot, [120,240], [47.5, 47.5], color=255, thick=2, linestyle=2
		oplot, [0,360], [mark_y, mark_y], color=255, thick=16, linestyle=0 
		oplot, [0,360], [mark_y, mark_y], color=0, thick=16, linestyle=2
	endif
	;ymin	= 0.8975362878E+00
	;oplot, [0, 360], [ymin, ymin], color=255, thick=16

	; ALTII-corners
	;oplot, [-22.5,157.5], [360-332.24, 360-332.24], color=255, thick=8, linestyle=2
	;oplot, [-22.5,157.5], [360-297.76, 360-297.76], color=255, thick=8, linestyle=2
	;oplot, [-22.5,157.5], [332.24, 332.24], color=255, thick=8, linestyle=2
	;oplot, [-22.5,157.5], [297.76, 297.76], color=255, thick=8, linestyle=2

	if (keyword_set(poincare_file)) then plot_poincare, poincare_file
	if (keyword_set(limiter_file))  then plot_limiter,  limiter_file, fill_color=fill_color, lim_style=lim_style
	if (keyword_set(grid_fileX))	then plot_grid, grid_fileX
	if (keyword_set(label))		then xyouts, lpos(0), lpos(1), label, charsize=1.5, charthick=2.0

	;xyouts, 315, 1.04, 'ALT-II target', charsize=1.5, charthick=2.0, orientation=90
	;xyouts, 320, 1.02, 'ALT-II', charsize=1.2, charthick=2.0, orientation=90


	; replot coordinate system if xrange and yrange are specified
	if (keyword_set(xrange) and keyword_set(yrange)) then $
	plot, xrange, yrange, /nodata, xstyle=xstyle, ystyle=xstyle, xmargin=[XL_Margin,XR_Margin], $
		 charsize=char_size, charthick=char_thick, /noerase, $
		 xticks=xticks, xtickv=xtickv, xtickinterval=xtickinterval, $
		 yticks=yticks, ytickv=ytickv, ytickinterval=ytickinterval

	; plot color bar
	bar_dummy	= transpose(levels#[1,1])
	colBarStr	= z_title

	; units label
	if (keyword_set(ztitle)) then u_title = ztitle

        if (not keyword_set(no_cb)) then begin
	contour, bar_dummy, [0,1], levels, /fill, levels=levels, $
		 ystyle=1, xstyle=4, xmargin=[XL_Margin_CBar,XR_Margin_Cbar], /noerase, $
		 ytitle=u_title, c_colors=colors, ymargin=[XD_Margin_CBar, XU_Margin_CBar], $
		 charsize=char_size, charthick=char_thick, ytickname=zticks, ytickinterval=ztickinterval
	endif


	close_device, ps_plot=ps_plot
end; plot_lc
;-----------------------------------------------------------------------
