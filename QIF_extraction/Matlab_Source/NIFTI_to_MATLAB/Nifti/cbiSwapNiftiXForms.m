function h=cbiSwapNiftiXForms(filename,newfilename,resetsform)
% h=cbiSwapNiftiXForms(filename [,newfilename, resetsform])
% 
% Swaps the qform and sform info of a file and saves the result.
% If a newfilename is specified, the result is saved to this file.
% If resetsform is nonzero, sets sfomr_code to 0 and sform44 to []

if (nargin<1)
    help(mfilename)
    return
end

[d,h]=cbiReadNifti(filename);


qf=h.qform44;
qfc=h.qform_code;
sf=h.sform44;
sfc=h.sform_code;

h.qform44=sf;
h.qform_code=sfc;

if (nargin<3 || ~resetsform)
    h.sform44=qf;
    h.sform_code=qfc;
else
    h.sform44=[];
    h.sform_code=0;
    h.srow_x(:)=0;
    h.srow_y(:)=0;
    h.srow_z(:)=0;
end


if (nargin>1)
    filename=newfilename;
end

cbiWriteNifti(filename,d,h)


