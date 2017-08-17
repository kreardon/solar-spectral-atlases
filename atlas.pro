;-------------------------------------------------------------
;+
; NAME:
;       READHEADER
; PURPOSE:
;       Read a FITS header.
; CATEGORY:
; CALLING SEQUENCE:
;       readheader,file,kwd,value,err=err
; INPUTS:
;       file  = FITS file name
; KEYWORD PARAMETERS:
;         ERR=err : error (= 0, OK) 
; OUTPUTS:
;       kwd   = keyword name list
;       value = values
; COMMON BLOCKS:
; NOTES:
; MODIFICATION HISTORY:
;       V. Andretta, 15 Dec, 1992
;-
;-------------------------------------------------------------

        pro readheader,file,kwd,value,err=err

	;---------  open file  -----------
        on_ioerror,ioerr
        get_lun,lun
        openr,lun,file
        record = assoc(lun,bytarr(80,36))

	;---------  read header  -----------
        kwd = strarr(1)
        value = strarr(1)
        recn = 0
        while max(strpos(kwd,'END     ')) do begin
          hdr = string(record(recn))
          kwd = [kwd,strmid(hdr,0,8)]
          value = [value,strmid(hdr,9,71)]
          recn = recn + 1
        endwhile
        kwd = kwd(1:*)
        value = value(1:*)
        err = abs(max(strpos(kwd,'SIMPLE  ')))

	;---------  close file and normal end  -----------
        close,lun
        free_lun,lun
        return

	;---------  end after error -----------
ioerr:  err = 1
        return
        end

;-------------------------------------------------------------
;+
; NAME:
;       KWDVAL
; PURPOSE:
;       Read a keyword value from a FITS header.
; CATEGORY:
; CALLING SEQUENCE:
;       v = kwdval(kwdlst,valuelst,kwd,type=type,err=err)
; INPUTS:
;       kwdlst   = keyword name list
;       valuelst = value list
;       kwd      = keyword name
; KEYWORD PARAMETERS:
;         TYPE=type : type of output : 
;                       1 : byte
;                       2 : integer
;                       3 : longword integer
;                       4 : floating point
;                       5 : double precision floating
;                       6 : complex floating
;                       7 : string
;         ERR=err : error (= 0, OK) 
; OUTPUTS:
;       v = returned value
; COMMON BLOCKS:
; NOTES:
; MODIFICATION HISTORY:
;       V. Andretta, 15 Dec, 1992
;-
;-------------------------------------------------------------

        function kwdval,kwdlst,valuelst,kwd,type=type,err=err

	;---------  check parameters and keywords  -----------
        err = n_params() lt 3  
        if not keyword_set(type) then type = 1
        err = err or (type lt 1) or (type gt 7)
        if err eq 1 then return,0
        kwdpos = $ 
          where(strtrim(strupcase(kwdlst)) eq strtrim(strupcase(kwd)),nkwd)
        if nkwd ne 0 then value = valuelst(kwdpos(0))
        case type of
          1: if nkwd ne 0 then kwdval=byte(value) else kwdval=byte(0)
          2: if nkwd ne 0 then kwdval=fix(value) else kwdval=0 
          3: if nkwd ne 0 then kwdval=long(value) else kwdval=0L
          4: if nkwd ne 0 then kwdval=float(value) else kwdval=0.
          5: if nkwd ne 0 then kwdval=double(value) else kwdval=0.D0
          6: if nkwd ne 0 then kwdval=complex(value) else kwdval=complex(0.)
          7: if nkwd eq 0 then kwdval='' $
                else begin
                  pos1 = strpos(value,"'")
                  pos2 = strpos(value,"'",pos1+1)
                  if ((pos1 ne -1) and (pos2 ne -1)) then begin
                    kwdval = strtrim(strmid(value,pos1+1,pos2-pos1-1))
                  endif else kwdval = ''
                endelse
        endcase 

	;---------  normal end  -----------
        return,kwdval
        end

;-------------------------------------------------------------
;+
; NAME:
;       GETHEADER
; PURPOSE:
;       Read essential data for the n-th atlas from its FITS header(s).
; CATEGORY:
; CALLING SEQUENCE:
;       getheader,n,err=err
; INPUTS:
;       n = atlas number
; KEYWORD PARAMETERS:
;         ERR=err : error (= 0, OK) 
; OUTPUTS:
; COMMON BLOCKS:
;       atlases : 
;         n_file     = number of files forming the atlas
;         file_atlas = name of the files forming the atlas
;       header :
;         nhrec   =  number of records forming the header
;         nwave   =  number of wavelength/wavenumber point in the file
;         wavemin =  lower wavelength/wavenumber
;         wavemax =  upper wavelength/wavenumber 
;         dwave   =  step
;         zero    =  offset for intensities/fluxes
;         scale   =  scaling factor
;         waveunit = wavelength/wavenumber unit name: Angstrom/cm-1
;         unit    =  intensity/flux unit name
; NOTES:
;       - Uses the following procedures/functions: 
;           readheader
;           kwdval
; MODIFICATION HISTORY:
;       V. Andretta, 15 Dec, 1992
;-
;-------------------------------------------------------------

        pro getheader,n,err=err

	;---------  common block definitions -----------
        ;---------  (see calling procedure)  -----------
        common atlases, n_file,file_atlas
        common header, nhrec,nwave, $
                       wavemin,wavemax,dwave, $
                       zero,scalef,waveunit,unit, $
                       author,institution

        for j = 0,n_file(n)-1 do begin
          readheader,file_atlas(n,j),kwd,value,err=err
          if (err gt 0) then begin
            print,'>>> Error reading file ',file_atlas(n,j)
            return
          endif
          nhrec(j) = n_elements(kwd)/36
          nwave(j) = kwdval(kwd,value,'naxis1',type=3,err=err)
          wavemin(j) = kwdval(kwd,value,'crval1',type=4,err=err)
          dwave(j) = kwdval(kwd,value,'cdelt1',type=4,err=err)
          pix1 = kwdval(kwd,value,'crpix1',type=3,err=err)
          wavemax(j) = wavemin(j)+(nwave(j)-pix1)*dwave(j)
          wavemin(j) = wavemin(j)+(1-pix1)*dwave(j)
          waveunit = strlowcase(kwdval(kwd,value,'ctype1',type=7,err=err))
          if (strpos(waveunit,'lambda') ne -1) then begin
            waveunit = 'A'
            wavemin(j) = wavemin(j)*1.E10
            wavemax(j) = wavemax(j)*1.E10
            dwave(j) = dwave(j)*1.E10
          endif else begin
            waveunit = 'cm-1'
            wavemin(j) = wavemin(j)*1.E-2
            wavemax(j) = wavemax(j)*1.E-2
            dwave(j) = dwave(j)*1.E-2
          endelse
          scalef(j) = kwdval(kwd,value,'bscale',type=4,err=err)
	  if err eq 1 then scalef(j) = 1.
          zero(j) = kwdval(kwd,value,'bzero',type=4,err=err)
	  if err eq 1 then zero(j) = 0.
          unit = strlowcase(kwdval(kwd,value,'bunit',type=7,err=err))
          author = kwdval(kwd,value,'observer',type=7,err=err)
          if err eq 1 then author = ''
          institution = kwdval(kwd,value,'obsvtory',type=7,err=err)
          if err eq 1 then institution = ''
        endfor

	;---------  end  -----------
        return
        end

;-------------------------------------------------------------
;+
; NAME:
;       GETDATA
; PURPOSE:
;       Read data in the n-th atlas from its FITS file(s).
; CATEGORY:
; CALLING SEQUENCE:
;       getdata,n,i1,i2,data,err=err
; INPUTS:
;       n = atlas number
;       i1 = position in the file(s) of the first data point
;       i2 = position in the file(s) of the last data point
; KEYWORD PARAMETERS:
;         ERR=err : error (= 0, OK) 
; OUTPUTS:
;       data = data values
; COMMON BLOCKS:
;       atlases : 
;         n_file     = number of files forming the atlas
;         file_atlas = name of the files forming the atlas
;       header :
;         nhrec   =  number of records forming the header
;         nwave   =  number of wavelength/wavenumber point in the file
;         wavemin =  lower wavelength/wavenumber
;         wavemax =  upper wavelength/wavenumber 
;         dwave   =  step
;         zero    =  offset for intensities/fluxes
;         scale   =  scaling factor
;         waveunit = wavelength/wavenumber unit name: Angstrom/cm-1
;         unit    =  intensity/flux unit name
; NOTES:
; MODIFICATION HISTORY:
;       V. Andretta, 15 Dec, 1992
;-
;-------------------------------------------------------------

        pro getdata,n,i1,i2,data,err=err

	;---------  common block definitions -----------
        ;---------  (see calling procedure)  -----------
        common atlases, n_file,file_atlas
        common header, nhrec,nwave, $
                       wavemin,wavemax,dwave, $
                       zero,scalef,waveunit,unit, $
                       author,institution

        data = intarr(1)

	;---------  go through the file(s) forming the atlas  -----------
        on_ioerror,ioerr
        for j = 0,n_file(n)-1 do begin

	;---------  open unit -----------
          err = 0
          get_lun,lun
          openr,lun,file_atlas(n,j)

	;---------  read only data in the range of interest  -----------
          if (i1(j) le i2(j)) then begin
            if (i2(j)-i1(j)) ge 200000L then $
              recordsize = nwave(j) $ ; longwords!
            else $
              recordsize = 1440  ; longwords!
            filedata = assoc(lun,intarr(recordsize),nhrec(j)*2880)
            rec1 = i1(j)/recordsize
            rec2 = i2(j)/recordsize
            tempdata = intarr(1)
            for k = rec1,rec2 do tempdata = [tempdata,filedata(k)]
            tempdata = tempdata(1:*)
            indx1 = i1(j) - rec1*recordsize
            indx2 = i2(j) - rec1*recordsize
            data = [data,tempdata(indx1:indx2)]
          endif

	;---------  close unit and go ahead  -----------
          close,lun
          free_lun,lun
        endfor
        data = data(1:*)

	;---------  swap (if necessary) and scale data  -----------
        byteorder,data,/ntohs
        data = float(data)
        k1 = 0
        for j = 0,n_file(n)-1 do if (i1(j) le i2(j)) then begin
          k2 = k1 + i2(j)-i1(j)
          data(k1:k2) = scalef(j)*data(k1:k2) + zero(j) 
          k1 = k1 + k2 + 1
        endif

	;---------  normal end  -----------
        return

	;---------  end after error -----------
ioerr:  err = 1
        return
        end

;-------------------------------------------------------------
;+
; NAME:
;       ATLAS
; PURPOSE:
;       Return a (solar) spectrum in a given spectral range, from 
;         atlases stored in FITS file.
; CATEGORY:
; CALLING SEQUENCE:
;       atlas, wave1, wave2, wave, spectrum, 
;              atlas_n=atlas_n, wavenumber=wvnum, list=lst, help=hlp
; INPUTS:
;       wave1 = lower limit of the spectral range
;       wave2 = upper limit of the spectral range
; KEYWORD PARAMETERS:
;         ATLAS_N=n   : atlas number
;         /WAVENUMBER : use wavenumbers (cm-1)
;         /LIST       : list available atlases. 
;	  /HELP       : show a parameter/keyword list.
; OUTPUTS:
;       wave     = wavelength/wavenumber points in the range wave1-wave2 
;       spectrum = (solar) spectrum in the range wave1-wave2 
; COMMON BLOCKS:
;       atlases : 
;         n_file     = number of files forming the atlas
;         file_atlas = name of the files forming the atlas
;       header :
;         nhrec   =  number of records forming the header
;         nwave   =  number of wavelength/wavenumber point in the file
;         wavemin =  lower wavelength/wavenumber
;         wavemax =  upper wavelength/wavenumber 
;         dwave   =  step
;         zero    =  offset for intensities/fluxes
;         scale   =  scaling factor
;         waveunit = wavelength/wavenumber unit name: Angstrom/cm-1
;         unit    =  intensity/flux unit name
; NOTES:
;       - Uses the following procedures/functions: 
;           getheader 
;           getdata
;       - If BSCALE and/or BZERO appear in the header they are applied to data. 
; MODIFICATION HISTORY:
;       V. Andretta, 15 Dec, 1992
;-
;-------------------------------------------------------------
 
        pro atlas, wave1, wave2, wave, spectrum, $
            atlas_n=atlas_n, wavenumber=wvnum, list=lst, help=hlp

	;---------  common block definitions -----------
        common atlases, n_file,file_atlas
        common header, nhrec,nwave, $
                       wavemin,wavemax,dwave, $
                       zero,scalef,waveunit,unit, $
                       author,institution

	;---------  available atlases  -----------
        n_atlas = 4                  ; number of atlases
        max_subfile = 2              ; max num. of files per atlas
        n_file = intarr(n_atlas)     ; num. of files of the atlas 
        file_atlas = strarr(n_atlas,max_subfile) ; name(s) of the file(s)

        atlas_location = File_Dirname(Routine_Filepath(/Either),/Mark)
        
        n_file(0) = 1
        file_atlas(0,0) = atlas_location + 'jun.fits'
        
        n_file(1) = 1
        file_atlas(1,0) = atlas_location + 'beckers.fits'
;
        n_file(2) = 2
        file_atlas(2,0) = atlas_location + 'kurucz_3a.fits'
        file_atlas(2,1) = atlas_location + 'kurucz_3b.fits'
;
        n_file(3) = 1                                  ; (IR int.)
        file_atlas(3,0) = atlas_location + 'kittpeak.fits'

	;---------  other common block variables  -----------
        nhrec = intarr(max_subfile)    ; number of heading records in the file
        nwave = lonarr(max_subfile)    ; number of points in the file
        wavemin = fltarr(max_subfile)  ; lower wavelength/wavenumber
        wavemax = fltarr(max_subfile)  ; upper wavelength/wavenumber 
        dwave = fltarr(max_subfile)    ; step
        zero  = fltarr(max_subfile)    ; offset for intensities/fluxes
        scalef = fltarr(max_subfile)    ; scaling factor
        waveunit = ' '                 ; unit name: Angstrom/cm-1
        unit = ' '                     ; int./flux unit name


	;---------  check parameters and keywords  -----------
        OK = 1
        OK = OK and (n_params() eq 4)
        OK = OK and (n_elements(wave1) eq 1) 
        OK = OK and (n_elements(wave2) eq 1) 
        OK = OK or keyword_set(lst)
        if (not keyword_set(atlas_n)) then atlas_n=0 else atlas_n=atlas_n-1
        if ((atlas_n lt 0) or (atlas_n ge n_atlas)) then begin
          print,'>>> Error: wrong atlas number'
          OK = 0
        endif

	;---------  help  -----------
        if ((not OK) or keyword_set(hlp)) then begin
	  print,' Return a (solar) spectrum stored in FITS file(s).'
          print,' Calling sequence:'
          print,'atlas, wave1, wave2, wave, spectrum'
          print,' Parameters:'
          print,'wave1  = lower limit of the spectral range.      in'
          print,'wave2  = upper limit of the spectral range.      in'
          print,'wave     = wavelength/wavenumbers.               out' 
          print,'spectrum = spectrum in the range wave1-wave2.    out' 
	  print,' Keywords:'
          print,'   ATLAS_N=n   : atlas number.' 
          print,'   /WAVENUMBER : use wavenumbers (cm-1)'
          print,'   /LIST       : list available atlases.'
	  print,' Notes:'
          print,'   If BSCALE and/or BZERO appear in the header,', $ 
                ' they are applied to data. '
	  return
	endif

	;---------  list available atlases  -----------
	if keyword_set(lst) then begin
          print,' Available atlases:'
          for i = 0,n_atlas-1 do begin
            getheader,i,err=err          ; read common block variables 
            if err gt 0 then begin
              print,'>>> Error reading the header of the atlas n.',i
              return
            endif
            frmt = '(A,I2,3A)'
            if (i eq 0) then defstr = ' (default)' else defstr = '' 
            unitname = strlowcase(strtrim(unit))
            if unitname eq 'relint' then unitname = 'relative intensity' $
              else if unitname eq 'relflx' then unitname = 'relative flux' $
                else unitname = 'unit = '+unitname
            print,format=frmt,'n.',i+1,defstr,', ',unitname
            frmt = '(2X,A)'
            descrpt = 'author(s): '
            if strlen(author) gt 0 then descrpt = descrpt+author  
            if strlen(institution) gt 0 then $
              descrpt = descrpt+' ('+institution+')'
            print,format=frmt,descrpt
            for j = 0,n_file(i)-1 do begin
              frmt = '(2X,A,F8.2,A,F8.2,3A,F5.3,2A)'
              print,format=frmt, $
                    'range = ',wavemin(j),'-',wavemax(j),' ',waveunit, $
                    ', step = ',dwave(j),' ',waveunit
            endfor
          endfor
          return
        endif

	;---------  read header parameters  -----------
        getheader,atlas_n,err=err
        if err gt 0 then begin
          print,'>>> Error reading the header of the atlas n.',atlas_n
          return
        endif

	;---------  check whether user wants same units  -----------
	;---------  of the atlas                         -----------
        wave1 = float(wave1)
        wave2 = float(wave2)
        inversion = keyword_set(wvnum) xor (waveunit ne 'A')
        if inversion then begin ; otherwise converts Angstroms to cm-1
          wave1 = 1.E8/wave1    ; (or viceversa, if the atlas is in cm-1
          wave2 = 1.E8/wave2    ;  and user wants Angstroms)
        endif

	;---------  exchanges extrema, if necessary  -----------
        tmpmin = min([wave1,wave2])
        tmpmax = max([wave1,wave2])

	;---------  trim range, if necessary  -----------
        tmpmin = max([tmpmin,min(wavemin(0:(n_file(atlas_n)-1)))]) 
        tmpmin = min([tmpmin,max(wavemax(0:(n_file(atlas_n)-1)))]) 
        tmpmax = max([tmpmax,min(wavemin(0:(n_file(atlas_n)-1)))]) 
        tmpmax = min([tmpmax,max(wavemax(0:(n_file(atlas_n)-1)))]) 
        wave1 = tmpmin
        wave2 = tmpmax

	;---------  find first and last index of the range(s)  -----------
        index1 = lonarr(n_file(atlas_n))
        index2 = lonarr(n_file(atlas_n))
        for j = 0,n_file(atlas_n)-1 do begin
          index1(j) = long((wave1-wavemin(j))/dwave(j)+0.5)
          index1(j) = max([index1(j),0L])
          index2(j) = long((wave2-wavemin(j))/dwave(j)+0.5)
          index2(j) = min([index2(j),nwave(j)-1])
        endfor

	;---------  fill wavelength/wavenumber vector WAVE  -----------
        wave = 0
        for j = 0,n_file(atlas_n)-1 do if (index1(j) le index2(j)) then begin
          indexes = lindgen(index2(j)-index1(j)+1)+index1(j)
          wave = [wave,wavemin(j)+indexes*dwave(j)]
        endif
        wave = wave(1:*)
        if inversion then wave=1.E8/wave
        if n_elements(wave) lt 2 then $
          print,'>>> Warning: values outside allowed range'

	;---------  read data in the range wave1-wave2 -----------
        getdata,atlas_n,index1,index2,spectrum,err=err
        if err gt 0 then begin
          print,'>>> Error reading the data from the atlas n.',atlas_n
          return
        endif
        if inversion then begin
          wave = rotate(wave,2)
          spectrum = rotate(spectrum,2)
        endif

	;---------  end  -----------
        return
        end
