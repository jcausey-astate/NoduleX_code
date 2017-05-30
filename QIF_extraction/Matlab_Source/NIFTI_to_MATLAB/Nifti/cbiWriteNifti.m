function [byteswritten,hdr]=cbiWriteNifti(varargin);
% [byteswritten,hdr]=cbiWriteNifti(filename,data,hdr [,prec,subset,short_nan] ) 
% [byteswritten,hdr]=cbiWriteNifti(filename,niftiim [,prec,subset,short_nan] ) 
%  Uses user-defined header struct 
%  (for second option, niftiim should be a struct with a .data field and a .hdr field)
%
% [byteswritten,hdr]=cbiWriteNifti(filename,data [,prec, subset,short_nan]);
%  Creates a new header
% 
% The header is always passed through cbiCreateNiftiHeader to ensure 
% consistency with the data.
%
%  prec:      Precision (dataype) of output. Default (if no header is specified) is 'float32'. Use [] for default.
%             Should be a recognized Matlab (or nifti) data type string.
%  subset:    4x1 cell array describing image subset to save. 1-offset (Matlab-style).
%             Only the following options are supported:
%             - to save a single z-slice (e.g. 4):  subset={[],[],4,[]}
%             - to save a single/multiple volumes of a time series (e.g. volumes 2-9):  subset={[],[],[],[2 9]}
%  short_nan: NaN handling for signed short (int16) data. If 1, will treat save NaNs as -32768 (smallest 
%             representable number) reserving -32767..32767 for scaled data; otherwise will save NaN as 0.
%             Default is 1 (use -32768 as NaN).
  
  
if (nargin<2)
    help(mfilename)
    error('(cbiWriteNifti) not enough input arguments!')
elseif (nargin>6)
    help(mfilename)
end
    
fname=varargin{1};
prec=[];
subset=[];
short_nan=1;
firstoptindex=3;
if (isstruct(varargin{2}))
    % niftiim input
    if (~isfield(varargin{2},'data') || ~isfield(varargin{2},'hdr'))
        error('(cbiWriteNifti) niftiim struct must have a data and a hdr field!')
    end
    data=varargin{2}.data;
    hdr=varargin{2}.hdr;
else
    % data/hdr input
    data=varargin{2};
    if (nargin>2)
        if (isstruct(varargin{3}))
            hdr=varargin{3};
            firstoptindex=4;
        else
            hdr=[];
        end
    else
        hdr=[];
    end
end
if (nargin>=firstoptindex)
    prec=varargin{firstoptindex};
    if (nargin>=firstoptindex+1)
        subset=varargin{firstoptindex+1};
        if (nargin>=firstoptindex+2)
            short_nan=varargin{firstoptindex+2};
        end
    end
end

[pathstr,imgname,hdrname,ext,compressext,iscompressed,hdriscompressed]=cbiParseNiftiFilename(fname,0);

if (~exist('subset') || isempty(subset))
    subset={[],[],[],[]};
    use_subset=0;
else
    use_subset=1;
    if (iscompressed)
        error('(cbiWriteNifti) Writing to a subset of a file only supported for non-compressed files!')    
    end
end
if (~exist('short_nan') || isempty(short_nan))
  short_nan=1;
end
  

hdr.img_name=fullfile(pathstr,[imgname compressext]);

switch (ext)
  case '.nii'  
    hdr.single_file=1;
    hdr.hdr_name=hdr.img_name;
    % This hack is necessary as sprintf on Mac OS X is broken (refuses to print \0 character)
    hdr.magic=char(zeros(1,4));
    hdr.magic(1:3)=sprintf('%s','n+1');
    hdr.vox_offset=352;
  case {'.hdr','.img'}
    hdr.single_file=0;
    hdr.hdr_name=fullfile(pathstr,hdrname);
    % This hack is necessary as sprintf on Mac OS X is broken (refuses to print \0 character)
    hdr.magic=char(zeros(1,4));
    hdr.magic(1:3)=sprintf('%s','ni1');
    hdr.vox_offset=0; % important!
  otherwise
    error('Not a valid NIFTI-1 file name extension. Legal values are .nii, .hdr, .img, .nii.gz/Z/zip, img.gz/Z/zip');
end  
if (~isfield(hdr,'endian'))
    hdr.endian='native';
end

  
% $$$   switch (ext)
% $$$  case '.nii'  
% $$$   hdr.single_file=1;
% $$$   hdr.hdr_name=fname;
% $$$   hdr.img_name=fname;
% $$$   hdr.magic=sprintf('%s\0','n+1');
% $$$   hdr.vox_offset=352;
% $$$  case {'.hdr','.img'}
% $$$   hdr.single_file=0;
% $$$   hdr.hdr_name=fullfile(pathstr,[bname '.hdr']);
% $$$   hdr.img_name=fullfile(pathstr,[bname '.img']);
% $$$   hdr.magic=sprintf('%s\0','ni1');
% $$$   hdr.vox_offset=0; % important!
% $$$  case '.gz','.Z' % zipped
% $$$   error('No support for zipped NIFTI-1 format under Matlab.');
% $$$  otherwise
% $$$   error('Not a valid NIFTI-1 file name extension. Legal values are .nii, .hdr, .img');
% $$$ end

% Ensure header matches data
if (exist('prec') & ~isempty(prec))
  hdr=cbiCreateNiftiHeader(hdr,data,'matlab_datatype',prec);
else
  hdr=cbiCreateNiftiHeader(hdr,data);
end

if (~strcmp(class(data),hdr.matlab_datatype))
  disp(['Warning: scaling data from ' class(data) ' to ' hdr.matlab_datatype]);
  disp('To avoid this, cast data to desired format before calling cbiWriteNifti, e.g.')
  disp('cbiWriteNifti(''myfilename'',int16(data),hdr,''int16'')')
end

% get hdr scaling factor and convert data if necessary
[data,hdr]=convertData(data,hdr,short_nan);  

% Write header
no_overwrite=0;
[hdr,zid]=cbiWriteNiftiHeader(hdr,fname,no_overwrite,hdr.single_file);

% Prepare to write data
if (~hdr.single_file)
  zid=zopen(hdr.img_name,'w',hdr.endian);
end

headerdim=hdr.dim(2:5); % Matlab 1-offset - hdr.dim(1) is actually hdr.dim(0)
headerdim(headerdim==0)=1; % Force null dimensions to be 1
is5D=0;
if (hdr.dim(6)>1)
  if (hdr.dim(5)>1)
    error('(cbiWriteNifti) No support for 5D data with multiple time points!');
  end
  is5D=1;
  headerdim(4)=hdr.dim(6);
  disp('(cbiWriteNifti) 5D data set detected');
end

loadSize=zeros(4,1);
for n=1:4
  if (length(subset{n})==1)
    subset{n}=[subset{n} subset{n}];    
  elseif (length(subset{n})>2)
    error('(cbiWriteNifti) subset should be a scalar or 2-vector');
  end
  if (isempty(subset{n}))
    loadSize(n)=headerdim(n);
    subset{n}=[1 loadSize(n)];
  else
    loadSize(n)=subset{n}(2)-subset{n}(1)+1;  
  end
end

if (any(loadSize>headerdim(1:4)))
  error('(cbiWriteNifti) subset index larger than image dimensions!');
elseif (any(loadSize(1:2)<headerdim(1:2)))
  error('(cbiWriteNifti) no support for saving subvolumes of data; only entire z-slices may be saved.');
end

%  Write emptiness so that we can move to the right offset. (This library does not currently support extensions, which would otherwise go here)
% 11/10/05 PJ
if (ztell(zid) < hdr.vox_offset)
  %  zwrite(zid,0,sprintf('integer*%d',(hdr.vox_offset-ftell(zid))));
  c=hdr.vox_offset-ztell(zid);
  if (zwrite(zid,zeros(c,1),'uint8')~=c)
      error('(cbiWriteNifti) error writing extension padding')    
  end
  if (ztell(zid) ~= hdr.vox_offset)
      error('(cbiWriteNifti) error writing extension padding')    
  end
end


% $$$ % Move to beginning of data (Not necessary as we have already written to corrct offset)
% $$$ status = zseek(zid,hdr.vox_offset,'bof');
% $$$ if status 
% $$$   ferror(zid)
% $$$ end
% $$$ 


dataSize=prod(loadSize);
writeFormat=hdr.matlab_datatype;
switch (hdr.matlab_datatype)    
 case 'binary'
  error('(cbiWriteNifti) No support for binary data')
 case 'complex64'
  writeFormat=float32;
 case 'complex128'
  writeFormat=float64;
 case 'RGB'
  error('(cbiWriteNifti) No support for RGB data');
 case {'complex256','float128'}
  error('(cbiWriteNifti) No support for 128-bit data on this platform!');
end

bytesPerElement=cbiSizeofNifti(writeFormat);
% For each time point:

% Move to correct location
readOrigin=sub2ind(headerdim(1:4)',subset{1}(1),subset{2}(1),subset{3}(1),subset{4}(1))-1; % now we're in C-land, hence 0-offset
											   
% Elements to write every time point
volSize=prod(headerdim(1:3));
readSize=prod(loadSize(1:3));
% Difference between volSize and readSize => offset to seek after every write
readOffset=volSize-readSize;
% Current position in data array
currPos=1;
if (strfind(hdr.matlab_datatype,'complex'))
  % Every voxel corresponds to two elements in the file
  readOrigin=readOrigin*2;
  readSize=readSize*2;
  readOffset=readOffset*2;
end
% Position file at first voxel
byteswritten=0;
if (use_subset)
    % Only supported for uncompressed files, hence we use fseek, not zseek (which only works for reads)
    fseek(zid,readOrigin*bytesPerElement,'cof');
end
for t=subset{4}(1):subset{4}(2)  
  % Extract current subset of data
  saveData=data(currPos:currPos+readSize-1);
  if (strfind(hdr.matlab_datatype,'complex'))
    % Separate complex data into real and imaginary parts
    realData=real(saveData);
    imagData=imag(saveData);
    saveData=zeros(2*prod(size(saveData)),1);
    saveData(1:2:readSize/2-1)=realData;
    saveData(2:2:readSize/2)=imagData;
  end

  % Write converted subset to file
  count=zwrite(zid,saveData,writeFormat);
  byteswritten=byteswritten+count;
  if (count~=readSize)
    zclose(zid);
    error(['(cbiWriteNifti) Error writing to file ' hdr.img_name]);
  end
  if (readOffset>0 & use_subset)
    % Only supported for uncompressed files, hence we use fseek, not zseek (which only works for reads)
    % fseek to next time point
    fseek(zid,readOffset*bytesPerElement,'cof');
  end
  currPos=currPos+readSize;
end

zclose(zid);

byteswritten=byteswritten*bytesPerElement;

return
  
function [data,hdr]=convertData(data,hdr,short_nan);
% Scales and shifts data (using hdr.scl_slope and hdr.scl_inter)
% and changes NaN's to 0 or MAXINT for non-floating point formats
% Returns hdr with scale factor changed (bug fixed 20060824)

% Calculate scale factor for non-floating point data
  switch (hdr.matlab_datatype)    
   case 'binary'
    error('unsupported format')
   case 'uint8'
    MAXINT=2^8-1;
   case {'uint16','ushort'}
    MAXINT=2^16-1;
   case {'uint32','uint'}
    MAXINT=2^32-1;
   case 'uint64'
    MAXINT=2^64-1;
   case 'int8'
    MAXINT=2^7-1;
   case {'int16','short'}
    MAXINT=2^15-1;
   case {'int32','int'}
    MAXINT=2^31-1;
   case 'int64'
    MAXINT=2^63-1;
   otherwise
    MAXINT=0;
  end
  if (MAXINT)
    hdr.scl_slope=max(data(:))/MAXINT;
  end
  
  % Scale and shift data if scale factor is nonzero
  if (~isnan(hdr.scl_slope) & hdr.scl_slope~=0)
    if (hdr.scl_slope~=1 | hdr.scl_inter~=0)
      data=double(data);
      data=(data-hdr.scl_inter)./hdr.scl_slope;
      switch (hdr.matlab_datatype)    
       case {'binary','uint8','uint16','short','ushort','uint32','uint','int','uint64','int8','int16','int32','int64'}
	data=round(data);
       otherwise
	% nothing
      end
    end
  end
    
  % Change NaNs for non-floating point datatypes
  switch (hdr.matlab_datatype)    
   case {'binary','uint8','uint16','ushort','uint32','uint','int','uint64','int8','int16','int32','int64'}
    data(isnan(data))=0;
   case {'int16','short'}
    if (short_nan)
      data(isnan(data))=-32768;
    else
      data(isnan(data))=0;
    end
  end

