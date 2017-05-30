function [pathstr,imgname,hdrname,ext,compext,iscomp,ishdrcomp]=cbiParseNiftiFilename(fname,checkfileexists)
% [pathstr,imgname,hdrname,ext,compext,iscomp,ishdrcomp]=cbiParseNiftiFilename(fname,checkfileexists)
%
%

if (~exist('checkfileexists','var') | isempty(checkfileexists))
    checkfileexists=0;
end

[pathstr,bname,ext]=fileparts(fname);
compext=[];
iscomp=0;
if (strcmp(ext,'.gz') | strcmp(ext,'.Z') | strcmp(ext,'.zip'))
    compext=ext;
    iscomp=1;
    [p,bname,ext]=fileparts(bname);
end

switch (ext)
 case '.nii'  
  imgname=[bname ext];
  hdrname=[bname ext];
  ishdrcomp=iscomp;
 case {'.hdr','.img'}
  imgname=[bname '.img'];
  hdrname=[bname '.hdr'];
  ishdrcomp=0;
  if (checkfileexists)
      % check that both exist  
      if (~exist(fullfile(pathstr,[hdrname]),'file'))
          error('(cbiParseNiftiFilename) Header file not found!');
      end
      if (~exist(fullfile(pathstr,[imgname compext]),'file'))
          error('(cbiParseNiftiFilename) Image file not found!');
      end
  end
 otherwise
  error('Not a valid NIFTI-1 file name extension. Legal values are .nii, .hdr, .img, .nii.gz/Z/zip, .img.gz/Z/zip');
end

