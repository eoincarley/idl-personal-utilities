pro mreadfits_shm, fitsfiles, iindex, data, xllp, yllp, nxp, nyp, $
   debug=debug, only_tags=only_tags, $
   xll=xll, yll=yll, nx=nx , ny=ny, destroy=destroy, outsize=outsize, $
   reuse_index=reuse_index, noscale=noscale, no_unmap_shm=no_unmap_shm
;+
;   Name: mreadfits_shm
;
;   Purpose: mreadfits analog using file/memory map (ITT shmmap)
;
;   Input Parameters:
;      fitsfiles - vector of fitsfiles
;      xllp, yllp, nxp, nyp - optional subfield in pixels 
;                           - positional equivilent of Keywords {xll,yll,nx,ny}
;
;   Output Parameters:
;      index [,data] ; ssw standard header-structures [, image data] 
;
;   Keyword Parameters:
;
;      xll - optional Lower left X pixel corrd (can be vector 1:1 fitsfiles)
;      yll - optional Lower left Y pixel coord ( ditto )
;      nx - #pixels in X 
;      ny - #pixels in Y  
;      outsize - optional desired output size (->rebin)
;      reuse_index (switch) - if set -and- input index exists, dont redo
;      noscale - if set, ignore scaling params (bscale/bzero...)
;      no_unmap_shm - if set, do not remove existing shared memory segments
;
;   History:
;      18-jul-2008 - S.L.Freeland - quasi general take on an SDO/JSOC routine
;      28-jul-2008 - S.L.Freeland - add subfield via XLL, YLL, NX, NY 
;      16-may-2009 - S.L.Freeland - add OUTSIZE
;      31-jul-2009 - S.L.Freeeland - add /REUSE_INDEX
;       7-aug-2009 - S.L.Freeland - add positional param option for subfield
;       8-mar-2010 - S.L.Freeland - add /NOSCALE
;      25-mar-2010 - S.L.Freeland - 16 bit UINT fix (bzero)
;      30-mar-2010 - S.L.Freeland - tweak header buff size calculation (hmi L0)
;      13-apr-2010 - S.L.Freeland - assure naxisN scalar 
;       1-jun-2010 - S.L.Freeland - handle primary image in 1st extension case (jsoc->imcopy "as is" qq?)
;      14-jul-2010 - S.L.Freeland - clear any mreadfits_shm shared areas 
;                                   (in  case user did a CTRL-C on previous)
;      15-oct-2010 - S.L.Freeland - avoid BYTEORDER issue w/4GB (ITT?)
;                                   expand BSCALE options
;      18-apr-2011 - S.L.Freeland - fix scaling/typing issue (HMI Cont for example)
;                                   add history/version
;      31-may-2011 - S.L.Freeland - Version 1.2; for SUBFOV, recalculate crpixN
;
;   Restrictions:
;      yes... gazilliions (but it does work on STEREO, AIA, & HMI on Macs/Linux)
;      assume all FITS files have consistant:
;         NAXIS1&2, BITPIX, number of FITS header 2880b blocks
;      Cannot currently specify OUTSIZE and subfield
;      Note: this routine works on uncompressed data ; mreadfits_tilecomp is ~analogous 
;      tile compressed fits reader - for SDO (aia/hmi) use, 'read_sdo.pro' is a wrapper
;      which will call the approriate reader (e.g. mreadfits_shm vs mreadfits_tilecomp)
;        
;-
version=1.2 
dodata=n_params() gt 2
debug=keyword_set(debug)
nf=n_elements(fitsfiles)
fsize=file_size(fitsfiles(0))
headbuff=bytarr(((fsize/2880) < 10)*2880)  ; assumed max of 10 header blocks ; ? 
openr,lun,/get_lun,fitsfiles(0),/swap_if_little
readu,lun,headbuff
free_lun,lun
ess=last_nelem(where_pattern(headbuff,'END')) ; find end of header block
if ess eq -1  then begin 
   box_message,'Cannot find end of header block... bailing'
   return
endif
nh=ess/2880 + 1 ; (ess mod 2880 ne 0) 
headbuff=headbuff(0:(nh*2880)-1)
bitpix=(where_pattern(headbuff,'BITPIX',bpcnt))(0)
naxis1=(where_pattern(headbuff,'NAXIS1',nx1cnt))(0)
naxis2=(where_pattern(headbuff,'NAXIS2',nx2cnt))(0)
next=get_fits_nextend(fitsfiles(0)) ; assume homologous in this regard
if next gt 0 then begin 
   bitpix=last_nelem(where_pattern(headbuff,'BITPIX',bpcnt))
   naxis1=last_nelem(where_pattern(headbuff,'NAXIS1',nx1cnt))
   naxis2=last_nelem(where_pattern(headbuff,'NAXIS2',nx2cnt))
endif

bscale=where_pattern(headbuff,'BSCALE',bscnt)
bzero=where_pattern(headbuff,'BZERO',bzcnt)
bblank=where_pattern(headbuff,'BLANK',bbcnt)

scaleit=bscnt gt 0
zeroit=bzcnt gt 0 
blankit=bbcnt gt 0
nobz=bzcnt
if bpcnt * nx1cnt *nx2cnt eq 0 then begin 
   box_message,'Cannot find at least one important FITS tag..., bailing'
   return
endif

bparr=[8,16,32,64,-32,-64] ; FITS BITPIX possibilities
tparr=[1,2,3,14,4,5] ; corresponding IDL data TYPE

reads,string(headbuff(bitpix+10:bitpix+50)),bpx
reads,string(headbuff(naxis1+10:naxis1+50)),nx1
reads,string(headbuff(naxis2+10:naxis2+50)),nx2
if scaleit then reads,string(headbuff(bscale+10:bscale+50)),bsc
if zeroit then reads,string(headbuff(bzero+10:bzero+50)),bzr
if blankit then reads,string(headbuff(bblank+10:bblank+50)),blk


dtype=where(bpx eq bparr,bpcnt)
if bpcnt eq 0 then begin
   box_message,'Unknow BITPIX... bailing'
   return
endif

datatype=tparr(dtype)
htemp=headbuff
htemp=reform(headbuff,80,36,nh)
hdx=data_chk(htemp,/nx)
index=make_array(hdx,nf,/byte,/nozero)
dtemp=make_array(nx1,nx2,type=datatype,/nozero)
;
np=n_params()  ; check for FOV positional params
if np gt 3 then xll=xllp
if np gt 4 then yll=yllp
if np gt 5 then nx=nxp
if np gt 6 then ny=nyp

subfov=n_elements(xll) gt 0 and n_elements(yll) gt 0 ; 
rebin=keyword_set(outsize)

if dodata then begin 
   if subfov then begin
      if n_elements(nx) eq 0 then nx1=256< (nx1-xll-1) else nx1=nx
      if n_elements(ny) eq 0 then nx2=256< (nx2-yll-1) else nx2=ny
      if n_elements(xll) eq 1 then xll=replicate(xll,nf)
      if n_elements(yll) eq 1 then yll=replicate(yll,nf)
   endif else begin 
      case n_elements(outsize) of
         0:
         1: begin  
               nx1=outsize(0)  
               nx2=outsize(0)
         endcase
         else:  begin 
            nx1=outsize(0) 
            nx2=outsize(1)
         endcase
      endcase
   endelse
   data=make_array(nx1,nx2,nf, type=datatype,/nozero)
   template={header:htemp,data:dtemp}
endif else template={header:htemp}
segs='mreadfits_shm'+strtrim(alphagen(nf),2)

help,/shared,out=inuse
ssmshm=where(strpos(inuse,'MREADFITS_SHM') ne -1,shmcnt)

if shmcnt gt 0 then begin 
   for s=0,shmcnt-1 do begin
      iuseg=ssw_strsplit(inuse(ssmshm(s)),' ',/head)
      box_message,'unmapping existing segment> '+ iuseg
      shmunmap,iuseg(0)
   endfor
endif


noindex=keyword_set(resuse_index) and data_chk(iindex,/struct) 
doindex=1-noindex
index=make_array(80,36,nh,nf,/byte)

for i=0,nf-1 do begin 
   shmmap,segs(i),get_name=gn,template=template,$
      filename=fitsfiles(i),/private
   delvarx,img
   img=shmvar(gn)
   header=img.header
   if dodata then begin
      ndata=temporary(img.data)
      case 1 of
         subfov: ndata=temporary(ndata(xll(i):xll(i)+nx1-1, yll(i):yll(i)+nx2-1))
         rebin: ndata=congrid(temporary(ndata),nx1,nx2)
         else:
      endcase
      data(0,0,i)=temporary(ndata)
   endif
   index(0,0,0,i)=header
   shmunmap,gn                                               
   if debug then help,/shared_memory
endfor

; 13-oct-2010 - 4GB problem with byteorder... 

if dodata and is_lendian() then begin ; handle arch byteorder issues...
   b4g=get_nbytes(data) ge 2.^32
   if b4g then begin 
      box_message,'avoiding 4GB byteorder limit'
      for im=0,nf-1 do data(0,0,im)=swap_endian(data(*,*,im))
   endif else swap_endian_inplace,data
endif
; optional scaling/bzero handling
upscale=0
if scaleit or zeroit and dodata then begin
   case 1 of 
      keyword_set(noscale): ; ignore scaling params
      bzr eq 2.^15 and bsc eq 1:data=temporary(uint(temporary(data))-bzr)
      bzr eq 2.^31 and bsc eq 1:data=temporary(ulong(temporary(data))-bzr)
      bzr eq 0. and bsc eq 1: ; no scaling...
      bsc ne 0: begin 
         if size(bsc,/tname) ne 'DOUBLE' then data *= float(bsc) else $
            data *= bsc
         if bzr ne 0 then begin 
            if size(bzr,/tname) ne 'DOUBLE' then data += float(bzr) else $
               data += bzr
         endif
         upscale=1
      endcase
      else: box_message,'Not yet trained to apply this scaling: bzr,bsc' + arr2str([bzr,bsc])
   endcase
endif

; header string->structure (via mreadfits_header)
if doindex then begin 
   t0=reltime(/now)
   headers=temporary(reform(string(index),36*nh*nf,/overwrite))
   mreadfits_header,fitsfiles,iindex,all_header=headers, only_tags=only_tags
   if upscale then begin  
      if required_tags(iindex,'bzero,bscale') then begin 
         iindex.bzero=0
         iindex.bscale=1.
         if tag_exist(iindex,'BLANK') then iindex.blank=blk*bsc+bzr
      endif
   endif
   update_history,iindex,/caller,version=version
endif
;
; date_obs fixup required
if required_tags(index,'date_obs,date_d$obs') then $ 
   if n_elements(all_vals(index.date_obs)) eq 1 then index.date_obs=index.date_d$obs

case 1 of 
   dodata and keyword_set(outsize): mreadfits_fixup,iindex,data
   dodata and subfov: begin ; adjust FITS FOV descriptors
      iindex.crpix1=iindex.crpix1-(xll-1)
      iindex.crpix2=iindex.crpix2-(yll-1)
      iindex.naxis1=nxp
      iindex.naxis2=nyp
      if tag_exist(iindex,'xcen') then $
        iindex.xcen=comp_fits_cen(iindex.crpix1,iindex.cdelt1, $
           gt_tagval(iindex,/crval1,missing=0.0))
      if tag_exist(iindex,'ycen') then $
        iindex.ycen=comp_fits_cen(iindex.crpix2,iindex.cdelt2, $
           gt_tagval(iindex,/crval2,missing=0.0))


   endcase
   else:
endcase


if debug then stop,'before return'
return
end



