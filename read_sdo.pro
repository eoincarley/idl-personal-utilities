pro read_sdo, files, index, data, xllp, yllp , nxp, nyp, _extra=_extra, $
   nodata=nodata, fnames_uncomp=fnames_uncomp, mixed_comp=mixed_comp, $
   parent_out=parent_out,outdir=outdir, comp_header=comp_header, $
   time_tag=time_tag
;+
;   Name: read_sdo
;
;   Purpose: read sdo/jsoc export files, aia and hmi
;
;   Input Parameters:
;      files - list of one or more sdo FITS files, jsoc/rice compressed or not
;      llx,lly,nx,ny - optionally, desired sub field (in pixels)
;
;   Output Parameters:
;      index - the ssw/fits meta data (structure vector)
;      data - optionally, the data (2D or 3D, depending on input files)
;
;   Keyword Parameters:
;      only_tags - (via inherit) - optional desired subset of tags/params
;      parent_out -  optional parent for uncompressed tilecomp
;      outdir - synonym for PARENT_OUT
;      outsize - (via inherit) - optional down-sizing request
;      _extra - unspecified keywords -> mreadfits_<blah> via inheritance
;               (see doc headers for mreadfits_shm & mreadfits_tilecomp, since
;               this should auto-track keyword/option evolution of those routines) 
;      nodata - (switch) - if set, do headers only, even if 3rd param 
;                           is included in call
;      fnames_uncomp - optional output keyword (compressed input files only)
;      UNCOMP_DELETE - (via inherit) If set, (and Compressed read), delete UnCompressed after reading
;      use_shared_lib (via inherit -> mreadfits_tilecomp) - if set - use memory-only
;         method if shared object available for OS/ARCH - Keh-Cheng Chu & Marc Derosa
;      mixed_comp - switch - if set, allow mixed compressed+uncompress logic
;      noshell - switch - avoid shell during uncompression stage (inherit->mreadfits_tilecomp)
;                  
;   History:
;      26-apr-2010 - S.L.Freeland - wrapper for mreadfits_shm/mreadfits_tilecomp
;      23-jun-2010 - S.L.Freeland - assures header only read for n_params=2
;                                   (ONLY_TAGS -> mreadfits_header.pro)
;      22-jul-2010 - S.L.Freeland - remove call_procedure - explicit comp/uncomp bifurcation
;      10-aug-2010 - S.L.Freeland - explicitly mention /UNCOMP_DELETE (-> mreadfits_tilecomp)
;       3-jan-2011 - S.L.Freeland - add /MIXED_COMP keyword+logic
;                    explicit PARENT_OUT + OUTDIR (replace inherit)
;      24-jan-2011 - S.L.Freeland - add NOSHELL document (inherit->mreadfits_tilecomp)
;      22-mar-2011 - S.L.Freeland - for compressed, return 
;                    "expected" values for some tags (bitpix,naxis1,naxis2...)
;                    (uncompressed header vals -> output index)
;                    Override with /COMP_HEADER switch
;       7-nov-2011 - S.L.Freeland - add TIME_TAG keyword (for orphaned jsoc files)
;       5-nov-2012 - S.L.Freeland - add /USE_SHARED_LIB blurb
;
;   NOTES - SUGGESTED you try either /NOSHELL -or- /USE_SHARED_LIB keywords
;           either should provide substantial speedup.
;           Removed full-disk restriction on /NOSHELL (eg., subfields OK)
;           If you try /NOSHELL and it breaks, 
;           please notify me: freeland@lmsal.com
; 
;   Restrictions:
;      as of today, cannot mix & match compressed and non-compressed
;      NOTE: for tile-compressed, currently writing intermediate decompressed versions
;      as of today, /USE_SHARED_LIB only for Mac & Linux 64bit idl
;      (see mreadfits_tilecomp header for more details)
;-
;   
if not file_exist(files(0)) then begin
   box_message,'IDL> read_sdo,<filelist>,index [,data,llpx,llpy,nx,ny] [,/noshell] [/use_shared]
   return
endif

noshell=keyword_set(noshell)

if keyword_set(mixed_comp) then begin 
   nfiles=n_elements(files)
   fsize=file_size(files)
   aac=where(fsize lt 33569280,ccnt)
   if ccnt gt 0 and ccnt ne nfiles then begin 
      box_message,'Mixed compression'
      ifiles=files
      read_sdo,files(aac),iizz,ddzz,/only_uncompress,fnames_uncomp=unames, $
          parent_out=parent_out, outdir=outdir, noshell=noshell
      files(aac)=unames
   endif ; else box_message,'/MIXED_COMP set but already homogenous'
endif
next=get_fits_nextend(files(0))

proc=(['mreadfits_shm','mreadfits_tilecomp'])(next gt 0)

nodata=keyword_set(nodata)

case 1 of 
   n_params() lt 2: box_message,'IDL> read_sdo,files,index [,data [,xll,yll,nx,ny]]
   n_params() eq 2 or nodata: mreadfits_header,files,index,exten=next,_extra=_extra
   n_params() eq 3: begin
      if next eq 0 then mreadfits_shm,files,index,data,_extra=_extra else $
                   mreadfits_tilecomp,files,index,data,_extra=_extra, fnames_uncomp=fnames_uncomp, $
                   parent_out=parent_out, outdir=outdir, time_tag=time_tag
   endcase
   else: begin
      if next eq 0 then mreadfits_shm,files,index,data,xllp,yllp,nxp,nyp,_extra=_extra else $
                   mreadfits_tilecomp,files,index,data,xllp,yllp,nxp,nyp,_extra=_extra, fnames_uncomp=fnames_uncomp, $
		      parent_out=parent_out, outdir=outdir
   endcase
endcase

comp2head=1-keyword_set(comp_header)
if next gt 0 and data_chk(index,/struct) and comp2head then begin 
   ftags=['BITPIX','NAXIS1','NAXIS2']
   ztags='Z'+ftags
   if required_tags(index,ftags) and required_tags(index,ztags) then begin 
      for i=0,n_elements(ftags)-1 do begin 
         index.(tag_index(index(0),ftags(i)))=gt_tagval(index,ztags(i))      
      end   
   endif
endif

if data_chk(index,/struct) and 1-tag_exist(index(0),'xcen') then begin 
   ; add xcen/ycen
   if required_tags(index(0),'crpix1,cdelt1') then begin 
      xcen=comp_fits_cen(index.crpix1,index.cdelt1,index.naxis1,index.crval1)
      ycen=comp_fits_cen(index.crpix2,index.cdelt2,index.naxis2,index.crval2)
      index=add_tag(index,xcen,'xcen')
      index=add_tag(index,ycen,'ycen')
   endif
endif

if n_elements(ifiles) eq n_elements(files) then begin
   if required_tags(_extra,'uncomp_delete') then begin 
      box_message,'removing mixed_comp uncompressed'
      ssw_file_delete,files(aac)
   endif
   files=ifiles ; restore input
endif

return
end


