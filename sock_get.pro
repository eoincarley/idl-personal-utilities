;+
; Project     : VSO
;
; Name        : SOCK_GET
;
; Purpose     : Wrapper around IDLnetURL object to download
;               files via HTTP and FTP
;
; Category    : utility system sockets
;
; Syntax      : IDL> sock_get,url,out_name,out_dir=out_dir
;
; Inputs      : URL = remote URL file name to download
;               OUT_NAME = optional output name for downloaded file
;
; Outputs     : LOCAL_FILE = full name of copied file
;
; Keywords    : OUT_DIR = output directory to download file
;               CLOBBER = clobber existing file
;               STATUS = 0/1/2 fail/success/file exists
;               CANCELLED = 1 if download cancelled
;               PROGRESS = show download progress
;
; History     : 27-Dec-2009, Zarro (ADNET) - Written
;                8-Oct-2010, Zarro (ADNET) - Dropped support for
;                COPY_FILE. Use LOCAl_FILE.
;              
;-

;-----------------------------------------------------------------  
function sock_get_callback, status, progress, data  

if (progress[0] eq 1) and (progress[1] gt 0) then begin
 if ptr_valid(data) then begin
  if (*data).progress then begin
   (*data).completed=progress[1] eq progress[2]
   if ~(*data).completed and ~(*data).cancelled then begin
    if ~widget_valid((*data).pid) then begin
     if allow_windows() then begin
      ourl=(*data).ourl
      if obj_valid(ourl) then begin
       ourl->getproperty,url_hostname=server,url_path=file
       bsize=progress[1]
       bmess=trim(str_format(bsize,"(i10)"))
       cmess=['Please wait. Downloading...','File: '+file,$
              'Size: '+bmess+' bytes',$
              'From: '+server,$
              'To: '+(*data).ofile]
        (*data).pid=progmeter(/init,button='Cancel',_extra=extra,input=cmess)
      endif
     endif
    endif
   endif
   if widget_valid((*data).pid) then begin
    val = float(progress[2])/float(progress[1])
    if val le 1 then begin
     if (progmeter((*data).pid,val) eq 'Cancel') then begin
      xkill,(*data).pid
      (*data).cancelled=1b
      return,0
     endif
    endif 
   endif

  endif
 endif
endif

;(*data).ourl->getproperty,verbose=verbose
;if verbose then print,status
   
if ptr_valid(data) then if ((*data).completed or (*data).cancelled) then xkill,(*data).pid
 
; return 1 to continue, return 0 to cancel  

return, 1
end

;-----------------------------------------------------------------------------
  
pro sock_get,url,out_name,clobber=clobber,local_file=local_file,$
  progress=progress,err=err,status=status,passive=passive,cancelled=cancelled,$
  out_dir=out_dir,username=username,password=password,_ref_extra=extra,verbose=verbose

err='' & status=0

if is_blank(url) then begin
 pr_syntax,'sock_get,url,out_dir=out_dir'
 return
endif

if ~since_version('6.4') then begin
 err='Requires IDL version 6.4 or greater.'
 message,err,/cont
 return
endif

verbose=keyword_set(verbose)
CATCH, errorStatus
IF (errorStatus NE 0) THEN BEGIN  
      CATCH, /CANCEL  
  
      ; Display the error msg in a dialog and at the IDL  
      ; command line.  
;      r = DIALOG_MESSAGE(!ERROR_STATE.msg, TITLE='URL Error', $  
;         /ERROR)  
      if verbose then PRINT, !ERROR_STATE.msg  
  
      ; Get the properties more details about the error and  
      ; display at the IDL command line.  
      ourl->GetProperty, RESPONSE_CODE=rspCode, $  
         RESPONSE_HEADER=rspHdr, RESPONSE_FILENAME=rspFn  
      if verbose then begin
       PRINT, 'rspCode = ', rspCode  
       PRINT, 'rspHdr= ', rspHdr  
       PRINT, 'rspFn= ', rspFn  
      endif
      goto,bail
 ENDIF  
  

cancelled=0b
local_file=''
clobber=keyword_set(clobber)
stc=parse_url(url)
file=file_break(stc.path)
path=file_break(stc.path,/path)+'/'
if is_blank(file) then begin
 err='File name not included in URL path.'
 message,err,/cont
 return
endif

;-- default copying file with same name to current directory

odir=curdir()
ofile=file
if is_string(out_name) then begin
 tdir=file_break(out_name,/path)
 if is_string(tdir) then odir=tdir 
 ofile=file_break(out_name)
endif
if is_string(out_dir) then odir=out_dir
if ~test_dir(odir,/verbose,err=err) then return
ofile=local_name(concat_dir(odir,ofile))
if ~clobber then begin
 if file_test(ofile) then begin
  message,'Local file '+ofile+' already exists. Use /clobber to recopy.',/cont
  local_file=ofile
  status=2
  return
 endif
endif

user='anonymous'
pass='nobody'
if is_string(stc.username) then user=stc.username
if is_string(stc.password) then pass=stc.password
if is_string(username) then user=username
if is_string(password) then pass=password
url_scheme=stc.scheme
if is_blank(url_scheme) then url_scheme='http'
active=0
if is_number(passive) then active=1-passive

;-- create pointer to pass to callback function 

callback_function=''
if keyword_set(progress) then callback_function='sock_get_callback'
ourl = obj_new('IDLnetUrl')
callback_data=ptr_new({ourl:ourl,ofile:ofile,pid:0l, progress:keyword_set(progress),cancelled:0b,completed:0b})
url_path=path+file
ourl->setproperty,callback_function=callback_function,url_scheme=url_scheme,$
                  callback_data=callback_data,verbose=verbose,$
                  url_host=stc.host,url_username=user,url_password=pass,$
                  url_path=url_path,_extra=extra,ftp_connection_mode=active

;-- start download 

result = oUrl->Get(file=ofile)  

;-- check what happened

bail: 

chk=file_test(ofile)
if is_string(result) then begin
 if chk then begin
  message,'Successfully downloaded - '+ofile,/cont 
  local_file=ofile
  status=1
 endif else begin
  err='Download failed.'
  message,err,/cont
  if chk then file_delete,ofile,/quiet
 endelse
endif else begin 
 if (*callback_data).cancelled then begin
  err='Download cancelled.' 
  cancelled=1b
 endif else $
  err=url+' not found on server.'
 message,err,/cont
 if chk then file_delete,ofile,/quiet
endelse

heap_free,callback_data

return & end  
