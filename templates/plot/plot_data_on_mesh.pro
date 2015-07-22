;- grid ----------------------------------------------------------------
pro iscrape, lun, iout
	str	= ''
	readf, lun, str
	print, strmid(str,2,80)
	iout	= long(strmid(str,32,10))
end; iscrape
pro rscrape, lun, rout
	str	= ''
	readf, lun, str
	print, strmid(str,2,80)
	rout	= double(strmid(str,32,10))
end; rscrape

pro read_grid, grid_file, xoffset=xoffset, yoffset=yoffset, swapxy=swapxy, units=units
	common	grid,	grid_id, mesh, n_x, n_y, x_title, y_title, $
			x_min, x_max, y_min, y_max

	; 1. initialize grid
	str	= ''
	openr, lun, grid_file, /get_lun
	readf, lun, str
	if (strmid(str,2,7) ne 'grid_id') then begin
		print, 'no grid type defined'
		stop
	endif
	grid_id	= fix(strmid(str,12,4))
	print, 'grid id = ', grid_id


	; 2. coordinate system, layout and fixed/dominant coordinate
	UNSTRUCTURED    = 1
	SEMI_STRUCTURED = 2
	STRUCTURED      = 3
	MESH_2D		= 4
	LOCAL		= 0
	CARTESIAN	= 1
	CYLINDRICAL	= 2
	TORUS		= 4
	; 2.1. get coordinates and layout from grid_id
	coordinates	= grid_id / 100
	grid_id		= grid_id - 100*coordinates
	layout		= grid_id / 10
	fixed_coord	= grid_id - 10*layout

	; 2.2. setup coordinate indices
	if (fixed_coord eq 1) then begin
		icoord1 = 2
		icoord2 = 3
	end else if (fixed_coord eq 2) then begin
		icoord1 = 3
		icoord2 = 1
	end else if (fixed_coord eq 0  or  fixed_coord eq 3) then begin
		icoord1 = 1
		icoord2 = 2
        end else begin
		print, 'error: invalid id = ', fixed_coord, ' for fixed coordinate!'
		stop
	endelse


	; 2.3. read reference parameters
	c_label	= make_array(4, /string)
	c_label(0)	= 'Plotting coordinate [a.u.]'
	c_label(1)	= ''
	c_label(2)	= ''
	c_label(3)	= ''
	case coordinates of
	CARTESIAN: begin
		c_label(1)	= 'x [cm]'
		c_label(2)	= 'y [cm]'
		c_label(3)	= 'z [cm]'
	end
	CYLINDRICAL: begin
		c_label(1)	= 'Radial Position [cm]'
		c_label(2)	= 'Vertical Position [cm]'
		c_label(3)	= 'Toroidal Angle [deg]'
	end
	TORUS: begin
		c_label(1)	= 'Minor Radius [cm]'
		c_label(2)	= 'Poloidal Angle [deg]'
		c_label(3)	= 'Toroidal Angle [deg]'
		rscrape, lun, R0
	end
	LOCAL: begin
	end
	endcase


	; 3. read grid resolution and nodes
	if (layout eq MESH_2D  and  fixed_coord > 0) then begin
		print, 'loading unstructured grid ...'
		print, '... using 1st and 2nd coordinate for plotting'
		iscrape, lun, n_x
		iscrape, lun, n_y
		rscrape, lun, z0

		mesh	= make_array(n_x, n_y, 2, /fl)
		tmp	= make_array(2, /fl)
		for j=0L,n_y-1 do begin
		for i=0L,n_x-1 do begin
			readf, lun, tmp
			mesh(i,j,*)	= tmp(*)
		endfor
		endfor
		x_title	= c_label(1)
		y_title	= c_label(2)

	end else begin
		print, 'error: invalid grid layout ', layout, '!'
		print, 'fixed_coord = ', fixed_coord
		stop
	endelse


	; 4. close grid file and post processing
	free_lun, lun



	; 5. units
	if (keyword_set(units)) then begin
		REAL_SPACE		= 1
		NORMALIZED_SPACE	= 2
		INDEX_SPACE		= 3
		if (units eq NORMALIZED_SPACE  or  units eq INDEX_SPACE) then begin
			wx	= 1.0
			wy	= 1.0
			if (units eq NORMALIZED_SPACE) then begin
				wx = 1.0 / (n_x - 1.0)
				wy = 1.0 / (n_y - 1.0)
			endif

			for j=0L,n_y-1 do begin
			for i=0L,n_x-1 do begin
				mesh(i,j,0)	= wx * i
				mesh(i,j,1)	= wy * j
			endfor
			endfor
		endif
	endif


	; 6.1. offset in x- and y- coordinates for plotting
        if (keyword_set(xoffset)) then begin
		mesh(*,*,0)	= mesh(*,*,0) + xoffset
	endif
        if (keyword_set(yoffset)) then begin
		mesh(*,*,1)	= mesh(*,*,1) + xoffset
	endif

	; 6.2. reverse plotting coordinates
	if (keyword_set(swapxy)) then begin
		stmp	= x_title
		x_title	= y_title
		y_title	= stmp

		rtmp		= mesh(*,*,0)
		mesh(*,*,0)	= mesh(*,*,1)
		mesh(*,*,1)	= rtmp
	endif

	; 6.3. upper and lower boundaries
	x_min	= min(mesh(*,*,0))
	x_max	= max(mesh(*,*,0))
	y_min	= min(mesh(*,*,1))
	y_max	= max(mesh(*,*,1))
end; read_grid
;-----------------------------------------------------------------------



;-----------------------------------------------------------------------
pro read_data, data_file
	common	data,	z_data, z_title, u_title
	common	grid

	result	= QUERY_ASCII(data_file, info)
	if (result ne 1) then begin
		print, 'error in QUERY_ASCII for file: ', data_file
		stop
	endif

	; number of rows / data lines
	n_data	= info.lines
	n_mesh	= (n_x-1)*(n_y-1)
	if (n_data ne n_mesh) then begin
		print, 'error: unexpected resolution!'
		print, 'mesh has ',n_x-1,' x ',n_y-1,' = ',n_mesh,' cells'
		print, 'data file has ',n_data,', rows'
		stop
	endif

	; number of columns
	columns	= info.words / n_data


	z_data	= make_array(n_x-1, n_y-1, columns, /fl)
	z_title	= ''
	u_title	= ''
	tmp	= make_array(columns, /fl)
	pi2	= 2*!PI
	openr, lun, data_file, /get_lun
	i	= 0L
	j	= 0L
	k	= 0L
	str	= ''
	while not eof(lun) do begin
		readf, lun, str
		if (strmid(str,0,1) ne '#') then begin
			reads, str, tmp
			z_data(i,j,*)	= tmp(*)
			i = i + 1
			k = k + 1
			if (i eq n_x-1) then begin
				i	= 0
				j	= j+1
			endif
		endif
	endwhile
	if (k ne n_data) then begin
		print, 'error: only ',k,' data lines found!'
		stop
	endif
	free_lun, lun
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



;-----------------------------------------------------------------------
pro plot_data_on_mesh, grid_file, data_file, idata, $
	xoffset=xoffset, yoffset=yoffset, swapxy=swapxy, $
	ps_plot=ps_plot, ct=ct, white0=white0, xsize=xsize, ysize=ysize, $
	xcells=xcells, ycells=ycells, $
	xrange=xrange, yrange=yrange, zrange=zrange, $
	units=units
	common	grid
	common	data
	common	cb_plot


	if (idata le 0) then begin
		print, 'error: invalid data column ', idata
		stop
	endif


	read_grid, grid_file, xoffset=xoffset, yoffset=yoffset, swapxy=swapxy, units=units


	read_data, data_file


	open_device, ps_plot=ps_plot, ct=ct, white0=white0, xsize=xsize, ysize=ysize


	if (keyword_set(zrange)) then begin
		z_min	= zrange(0)
		z_max	= zrange(1)
	endif else begin
		z_min	= min(z_data(*,*,idata-1))
		z_max	= max(z_data(*,*,idata-1))
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


	; x-range
	if (keyword_set(xrange)) then begin
		xrange_plot	= xrange
	endif else begin
		xrange_plot	= [x_min, x_max]
	endelse
	; y-range
	if (keyword_set(yrange)) then begin
		yrange_plot	= yrange
	endif else begin
		yrange_plot	= [y_min, y_max]
	endelse

	; initialize plot window
	plot, xrange_plot, yrange_plot, xstyle=1, ystyle=1, /nodata, $
		xmargin=[XL_Margin, XR_Margin]

	; cell range
	i1	= 0L
	i2	= n_x-2
	j1	= 0L
	j2	= n_y-2
	if (keyword_set(xcells)) then begin
		if (xcells(0) ge 0  and  xcells(0) le n_x-2) then begin
			i1	= xcells(0)
		endif else begin
			print, 'error: lower cell boundary out of range!'
			print, '0 <= xcells(0) <= ', n_x-2, ' required!'
			stop
		endelse
		if (xcells(1) ge 0  and  xcells(1) le n_x-2) then begin
			i2	= xcells(1)
		endif else begin
			print, 'error: upper cell boundary out of range!'
			print, '0 <= xcells(1) <= ', n_x-2, ' required!'
			stop
		endelse
	endif
	if (keyword_set(ycells)) then begin
		if (ycells(0) ge 0  and  ycells(0) le n_y-2) then begin
			j1	= ycells(0)
		endif else begin
			print, 'error: lower cell boundary out of range!'
			print, '0 <= ycells(0) <= ', n_y-2, ' required!'
			stop
		endelse
		if (ycells(1) ge 0  and  ycells(1) le n_y-2) then begin
			j2	= ycells(1)
		endif else begin
			print, 'error: upper cell boundary out of range!'
			print, '0 <= ycells(1) <= ', n_y-2, ' required!'
			stop
		endelse
	endif

	xarr	= make_array (4, /fl)
	yarr	= make_array (4, /fl)
	for i=i1,i2 do begin
	for j=j1,j2 do begin
		xarr(0)	= mesh(i,  j,  0)
		xarr(1)	= mesh(i+1,j,  0)
		xarr(2)	= mesh(i+1,j+1,0)
		xarr(3)	= mesh(i,  j+1,0)

		yarr(0)	= mesh(i,  j,  1)
		yarr(1)	= mesh(i+1,j,  1)
		yarr(2)	= mesh(i+1,j+1,1)
		yarr(3)	= mesh(i,  j+1,1)


		znn	= round((n_colors-1) * (z_data(i,j,idata-1) - z_min) / (z_max - z_min)) + 1
		if (z_data(i,j,idata-1) ge z_max) then znn = n_colors
		if (z_data(i,j,idata-1) le z_min) then znn = 1

		polyfill, xarr, yarr, col=znn, /data, noclip=0
	endfor
	endfor


	; color bar
	bar_dummy	= transpose(levels#[1,1])
	colBarStr	= z_title

	contour, bar_dummy, [0,1], levels, /fill, levels=levels, $
		 ystyle=1, xstyle=4, xmargin=[XL_Margin_CBar,XR_Margin_Cbar], /noerase, $
		 ytitle=utitle, c_colors=colors, $
		 charsize=char_size, charthick=char_thick, $
		 ytickname=zticks, ytickinterval=ztickinterval


	close_device, ps_plot=ps_plot
end; plot_data_on_mesh
;-----------------------------------------------------------------------
