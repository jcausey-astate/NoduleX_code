function varargout=cbiResliceNifti(inputfile,formatfile,outputfile,usesform)
% Reslices a Nifti file and adjusts qform/sform accordingly.
%
% 1. Reslice file to match another in same space (qform or sform):
%    [output] = cbiResliceNifti(inputfile,formatfile [,outputfile,usesform])
%  formatfile : Filename or Nifti header struct of target format file
%    usesform : Transform with respect to sform (default transforms wrt qform)
%
% 2. Extract a subvolume of a file
%    [output] = cbiResliceNifti(inputfile, subset, [,outputfile])
%      subset : cell array specifying subvolume to extract (3 or 4 long)
%               {[xmin xmax] [ymin ymax] [zmin zmax] [tmin tmx]}
%
% 3. Apply a 4x4 transformation mtx 
%    [output] = cbiResliceNifti(inputfile, mtx [,outputfile])
%         mtx : 4x4 homogeneous transformation matrix (world_coord_new=mtx*world_coord_old)
%    usesform : Transform with respect to sform (default transforms wrt qform)
%
% 4. Resample within same space
%    [output] = cbiResliceNifti(inputfile, res [,outputfile])
%         res : [i_res j_res k_res t_res] (3- or 4-vector)
%               Not that i,j,k refer to the array coordinates; these will only correspond to x,y,z for images
%               acquired in cardinal planes and ordered so that off-diagonal qform44 entries are zero.
%               Use 0 to indicate no change, eg [1 0 0 0] will resample to 1mm in i (x) 
%               and keep other dimensions
%               Negative values are interpreted as scalings, so [-2 -2 -2] means 
%               subsample 2x in each of i,j,k (x,y,z)
% 
%   inputfile : filename or Nifti data/hdr struct
%  outputfile : filename for output. Must be specified if no output is given.
%
%      output : Nifti data/hdr struct

%%    modify to work also on 4D+ data
%%    modify to reslice and rotate a subset of a volume (eg get a slice for overlay)
%%    modify to allow option to specify size, offset, rotation
%%    modify to allow option to specify explicit 4x4 transform

if (nargin<2)
    if (islinux)
        if (nargin<1)
            [status,inputfile]=system('zenity --title "Choose input file" --file-selection --file-filter="*.img *nii *img.gz *nii.gz" 2> /dev/null');
            if (isempty(inputfile)) return;end
        end
        [status,formatfile]=system('zenity --title "Choose format file" --file-selection --file-filter="*.img *nii *img.gz *nii.gz" 2> /dev/null')
        if (isempty(formatfile)) return;end
        inputfile=inputfile(1:end-1);
        formatfile=formatfile(1:end-1);
    else        
        help(mfilename)
        return
    end
end

if (~exist('usesform','var'))
    usesform=0;
end

% Load input file or use preloaded Nifti struct
if (ischar(inputfile))
    inp=cbiReadNifti(inputfile);
elseif (isstruct(inputfile))    
    inp=inputfile;
else
    help(mfilename)
    return
end

option=0;
% Load format (file, Nifti struct, or format cell array)
if (ischar(formatfile))
    fmt=cbiReadNifti(formatfile);
    option=1;
elseif (isstruct(formatfile))
    fmt=formatfile;
    option=1;
elseif (iscell(formatfile))
    % should check that it is a 3-long cell with ranges
    fmt=formatfile;
    if (length(fmt)<3 || length(fmt)>4 || any(cellfun('size',fmt,2)~=2))
        error('subset must be a cell array of 2 long row vectors')
    end
    option=2;
elseif (isnumeric(formatfile))
    fmt=formatfile;        
    if (all(size(formatfile)==[4 4]))
        option=3;
    elseif (isvector(fmt) && length(fmt)>= 3 &&  length(fmt)<=4)
        option=4;
    end        
end

if (option==1)
    xs=fmt.hdr.dim(2);
    ys=fmt.hdr.dim(3);
    zs=fmt.hdr.dim(4);
    formatdatasize=size(fmt.data);
    % Get array coordinates of resliced volume (format)
    %% NOTE: meshgrid and interp3 use Matlab's dumbass row major order, so x and y need to be switched in both function calls.
    %    [iy,ix,iz]=meshgrid(0:ys-1,0:xs-1,0:zs-1);
    [ix,iy,iz]=ndgrid(0:xs-1,0:ys-1,0:zs-1);

    % Compute transformation xform that maps from target array coordinates to input array coordinates
    % This is then used below to compute input array coordinates corresponding to the target array coords
    if (usesform)
        xform=inv(inp.hdr.sform44)*fmt.hdr.sform44;
    else
        xform=inv(inp.hdr.qform44)*fmt.hdr.qform44;
    end
elseif (option==2)
    % subset cell array describing volume to extract
    xrange=fmt{1}(1):fmt{1}(2);
    yrange=fmt{2}(1):fmt{2}(2);
    zrange=fmt{3}(1):fmt{3}(2);
    %    [iy,ix,iz]=meshgrid(yrange,xrange,zrange);
    [ix,iy,iz]=ndgrid(xrange,yrange,zrange);
    %formatdatasize=[diff(fmt{2}) diff(fmt{1}) diff(fmt{3})]+1;
    % xform is given by translation relative to original volume    
    coordshift=ones(4,1);
    coordshift(1)=fmt{1}(1)-1; % note subset is 1-offset
    coordshift(2)=fmt{2}(1)-1; 
    coordshift(3)=fmt{3}(1)-1;
    xform=[];
elseif (option==3)
    % compute new world coordinates
    mtx=fmt;
    % mtx*s/form is the new mapping from array to world coordinates (equivalent to fmt.hdr.q/sform in option 1)
    if (usesform)
        xform=inv(inp.hdr.sform44)*mtx*inp.hdr.sform44;
    else
        xform=inv(inp.hdr.qform44)*mtx*inp.hdr.qform44;
    end
    xs=inp.hdr.dim(2);
    ys=inp.hdr.dim(3);
    zs=inp.hdr.dim(4);
    %[iy,ix,iz]=meshgrid(0:ys-1,0:xs-1,0:zs-1);
    [ix,iy,iz]=ndgrid(0:xs-1,0:ys-1,0:zs-1);
elseif (option==4)
    % Determine xform
    oldres=inp.hdr.pixdim(2:4);
    resscale=fmt(1:3);
    resscale=resscale(:);
    % Convert absolute voxel dimensions to scaling parameters
    resscale(resscale>0)=resscale(resscale>0)./oldres(resscale>0);
    resscale(resscale==0)=1;
    resscale=abs(resscale);    
    mtx=diag([resscale;1]); 
    % Add translation parameters (difference between half voxels)
    newres=oldres.*resscale;
    coordshift=(newres-oldres)/2;
    xform=mtx;
    mtx(1:3,4)=coordshift;
    % Compute grid
    newsz=inp.hdr.dim(2:4)./resscale;
    xs=newsz(1);
    ys=newsz(2);
    zs=newsz(3);
    [ix,iy,iz]=ndgrid(0:xs-1,0:ys-1,0:zs-1);
end

fmtcoords=[ix(:) iy(:) iz(:) ones(prod(size(ix)),1)]';

% Transform
if (~isempty(xform))
    inpcoords=xform*fmtcoords;
else
    inpcoords=fmtcoords;
end

% Interpolate (note 1-offset is required here)
%% NOTE: meshgrid and interp3 use Matlab's dumbass row major order, so x and y need to be switched in both function calls.
if (isempty(xform))
    outp.data=inp.data(xrange,yrange,zrange);%reshape(interp3(inp.data,inpcoords(2,:),inpcoords(1,:),inpcoords(3,:),'nearest',NaN),size(ix));
else
    %    outp.data=reshape(interp3(inp.data,inpcoords(2,:)+1,inpcoords(1,:)+1,inpcoords(3,:)+1,'linear',NaN),size(ix));
    outp.data=reshape(interpn(inp.data,inpcoords(1,:)+1,inpcoords(2,:)+1,inpcoords(3,:)+1,'linear',NaN),size(ix));
end
outp.hdr=inp.hdr;
if (option==1)
    outp.hdr.qform44=fmt.hdr.qform44;
    outp.hdr.sform44=fmt.hdr.sform44;
    outp.hdr.qform_code=fmt.hdr.qform_code;
    outp.hdr.sform_code=fmt.hdr.sform_code;
    outp.hdr.dim=fmt.hdr.dim;
elseif (option==2)
    outp.hdr.qform44(:,4)=inp.hdr.qform44*coordshift;
    outp.hdr.sform44(:,4)=inp.hdr.sform44*coordshift;
    outp.hdr.qform44(4,4)=1;
    outp.hdr.sform44(4,4)=1;
elseif (option==3)
    % transform takes place in world coordinates, so second op
    outp.hdr.qform44=mtx*inp.hdr.qform44;
    outp.hdr.sform44=mtx*inp.hdr.sform44;    
elseif (option==4)
    % scaling takes place in array coordinates, so first op
    outp.hdr.qform44=inp.hdr.qform44*mtx;
    outp.hdr.sform44=inp.hdr.sform44*mtx;    
end

outp.hdr.dim(2:4)=size(outp.data);

if (nargout==1)
    varargout{1}=outp;
else
    if (nargin<3)
        [status,outputfile]=system('zenity --title "Choose output file" --file-selection --confirm-overwrite --file-filter="*.img *nii *img.gz *nii.gz" 2> /dev/null');
        if (isempty(outputfile)) return;end
        outputfile=outputfile(1:end-1);

    end
    cbiWriteNifti(outputfile,outp);
end

    