function varargout = cbiReadNiftiHeader(fname,tfi2nifti)
% hdr = cbiReadNiftiHeader(fname [,tfi2nifti])         Reads header and closes file
% [hdr,zid] = cbiReadNiftiHeader(fname [,tfi2nifti])   Reads header and returns file pointer for single files only
% 
% Loads the header from a NIFTI-1 file.  
% The header struct contains all stored fields plus 
% some additional extracted fields. These are:
% .hdr_name    - name of header file
% .img_name    - name of image file
% .is_analyze  - 1 if valid Analyze, but not a valid NIFTI file, 0 otherwise
% .single_file - flag for single file (1 for .nii, 0 for .hdr/.img)
% .endian      - endianness of file
% .matlab_datatype - data type in Matlab string format (see fread)
% .qform44    - homogeneous transform matrix extracted from qform data (rigid body - native coordinate system) ([] for Analyze) 
% .sform44    - homogeneous transform matrix extracted from sform data (affine - optional alignment transform) ([] for Analyze) 
%
% fname can have any legal extension (.hdr, .img, .nii, .img.gz, nii.gz).
%
% Optional second argument tfi2nifti (default tfi2nifti=1):
% - If nonzero, will convert TFI (SurfRelax) Analyze origin and aspect info to a Nifti-style qform44 and set quaternion parameters.
%   This allows using TFI Analyze files transparently as Nifti files (with the limitation that the TFI Analyze header only represents
%   scalings and translations, not rotations). Specifically, it allows running mrAlign with TFI Analyze files as destination files
%   and the derived xform can then be used in VolumeViewer.

% Check name and file type.
[pathstr,imgname,hdrname,ext,compressext,iscompressed,hdriscompressed]=cbiParseNiftiFilename(fname,1);

hdr.img_name=fullfile(pathstr,[imgname compressext]);

if (nargin<2)
    tfi2nifti=1;
end

switch (ext)
  case '.nii'  
    hdr.single_file=1;
    hdr.hdr_name=hdr.img_name;
  case {'.hdr','.img'}
    hdr.single_file=0;
    hdr.hdr_name=fullfile(pathstr,hdrname);
  otherwise
    error('Not a valid NIFTI-1 file name extension. Legal values are .nii, .hdr, .img');
end

% Find the endian-ness of the file
hdr.endian='b'; % open file in big-endian
zid=zopen(hdr.hdr_name,'r',hdr.endian);

% check if this gives the correct header size - if not use little-endian
testval = zread(zid,1,'int32');
if (testval ~= 348)
  zclose(zid);
  hdr.endian='l';
  zid=zopen(hdr.hdr_name,'r',hdr.endian);
  testval = zread(zid,1,'int32');
  if (testval ~= 348)
    error('Illegal file format - incorrect header size (should be 348 bytes)');
  end
end

%% --- was header_key substruct ---
hdr.sizeof_hdr = testval;
dummy = zread(zid,35,'char');	% ditch the remaining initial header stuff
hdr.dim_info = zread(zid,1,'char');
%% --- was image_dimension substruct ---
hdr.dim = zread(zid,8,'int16');
hdr.intent_ps = zread(zid,3,'float');
hdr.intent_code = zread(zid,1,'int16');
hdr.datatype = zread(zid,1,'int16');
hdr.bitpix = zread(zid,1,'int16');
hdr.slice_start = zread(zid,1,'int16');
hdr.pixdim = zread(zid,8,'float');
hdr.vox_offset = zread(zid,1,'float');
hdr.scl_slope = zread(zid,1,'float');
hdr.scl_inter = zread(zid,1,'float');
hdr.slice_end = zread(zid,1,'int16');
hdr.slice_code = zread(zid,1,'char');
hdr.xyzt_units = zread(zid,1,'char');
hdr.cal_max = zread(zid,1,'float');
hdr.cal_min = zread(zid,1,'float');
hdr.slice_duration = zread(zid,1,'float');
hdr.toffset = zread(zid,1,'float');
dummy = zread(zid,2,'int32');
%% --- was data_history substruct ---
hdr.descrip = char(transpose(zread(zid,80,'char')));
hdr.aux_file = char(transpose(zread(zid,24,'char')));

hdr.qform_code = zread(zid,1,'int16');
hdr.sform_code = zread(zid,1,'int16');
hdr.quatern_b = zread(zid,1,'float');
hdr.quatern_c = zread(zid,1,'float');
hdr.quatern_d = zread(zid,1,'float');
hdr.qoffset_x = zread(zid,1,'float');
hdr.qoffset_y = zread(zid,1,'float');
hdr.qoffset_z = zread(zid,1,'float');
hdr.srow_x    = zread(zid,4,'float')';
hdr.srow_y    = zread(zid,4,'float')';
hdr.srow_z    = zread(zid,4,'float')';
hdr.intent_name  = char(transpose(zread(zid,16,'char')));
hdr.magic        = char(transpose(zread(zid,4,'char')));

if (nargout<2 | ~hdr.single_file)
    zclose(zid);
end

% Set data type
hdr.matlab_datatype=cbiNiftiDatatype2Matlab(hdr.datatype);

hdr.is_analyze=(~hdr.single_file & isempty(strfind(hdr.magic,'ni1')));

% Extract xform and sform
if (hdr.is_analyze)
    if (tfi2nifti)
        xform=eye(4);
        xform(1:3,1:3)=diag(hdr.pixdim(2:4));
        xform(1:3,4)=-hdr.pixdim(6:8);
        hdr.sform44=[];
        hdr=cbiSetNiftiQform( hdr, xform );
    else
        hdr.qform44=[];
        hdr.sform44=[];
    end
else
  hdr.qform44=cbiQuaternionToHomogeneous(hdr);
  hdr.sform44=eye(4);
  if (hdr.sform_code>0)
    hdr.sform44(1,:)=hdr.srow_x;
    hdr.sform44(2,:)=hdr.srow_y;
    hdr.sform44(3,:)=hdr.srow_z;
  end
end

varargout{1}=hdr;
if (nargout==2)
    if (hdr.single_file)
        varargout{2}=zid;
    else
        varargout{2}=[];
    end
end
