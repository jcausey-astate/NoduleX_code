% cbiGetSwapVect.m
%
%      usage: cbiGetSwapVect(filename1,filename2)
%         by: justin gardner; modified by Jonas Larsson for Nifti implementation
%       date: 06/16/05
%    purpose: find which dimensions to swap or flip
%             so that the inplane anatomy and off-line
%             reconstructed epi images match
%             will return a swap vector which can then be passed
%             to cbiSwapNiftiDimensions to swap data and Nifti header
%              
%       e.g.: swapdim('Anatomy/Inplane/jg050614+03+t1_2D_mprage_24sl.hdr','Pfiles_preMC/09+cbi_bold2_1500_3x3_24sl_recon.hdr')
%
%             Note that unlike avwswapdim, cbiSwapNiftiDimensions 
%             correctly modifies header information to preserve 
%             rotation fields in qform and sform.
%             THEREFORE DO NOT USE avwswapdim. 
%
%             When called with no arguments, looks for
%             the anatomy file in Raw/Anatomy/Inplane
%             and the epi from either Raw/Pfiles_preMC or
%             Raw/Pfiles
%
%       e.g.: swapdim;
%
%
function [cbiswapvect, swapname] = cbiGetSwapVect(filename1,filename2)

avwcommand = [];swapname = [];
if (nargin == 0)
  % see if we can find the anatomy and an epi
  filename2 = getimgfromdir('Raw/Anatomy/Inplane');
  if isempty(filename2),help swapdim;return,end
  if (isdir('Raw/Pfiles_preMC'))
    filename1 = getimgfromdir('Raw/Pfiles_preMC');
  else
    filename1 = getimgfromdir('Raw/Pfiles');
  end
  if isempty(filename1),help swapdim;return,end
elseif (nargin == 1)
  % see if we can find the anatomy and an epi
  filename2 = getimgfromdir(sprintf('%s/Raw/Anatomy/Inplane',filename1));
  if isempty(filename2),help swapdim;return,end
  if (isdir(sprintf('%s/Raw/Pfiles_preMC',filename1)))
    filename1 = getimgfromdir(sprintf('%s/Raw/Pfiles_preMC',filename1));
  else
    filename1 = getimgfromdir(sprintf('%s/Raw/Pfiles',filename1));
  end
  if isempty(filename1),help swapdim;return,end
elseif (nargin ~= 2)
  help swapdim;
  return
end

swapname = stripext(filename1);

% if there is no file extension then add it
if (isempty(getext(getlastdir(filename1))))
  filename1 = sprintf('%s.hdr',filename1);
end  
if (isempty(getext(getlastdir(filename2))))
  filename2 = sprintf('%s.hdr',filename2);
end  

% check to see if we have the right file extension
if (~sum(strcmp(getext(getlastdir(filename1)),{'hdr' 'img'})))
  disp(sprintf('UHOH: %s does not have correct extension (.hdr or .img)',filename1));
  return
end
if (~sum(strcmp(getext(getlastdir(filename2)),{'hdr' 'img'})))
  disp(sprintf('UHOH: %s does not have correct extension (.hdr or .img)',filename2));
  return
end

% make sure file exists 
if (~isfile(filename1))
  disp(sprintf('UHOH: Could not find file %s',filename1));
  return
end
if (~isfile(filename2))
  disp(sprintf('UHOH: Could not find file %s',filename2));
  return
end

% open the headers
hdr1 = cbiReadNiftiHeader(filename1);
hdr2 = cbiReadNiftiHeader(filename2);

% get the rotation matrices
R1 = quaternion2rotation(hdr1);
R2 = quaternion2rotation(hdr2);

% and compute the transformation that will move R1 to R2
transform = (R1') * ((R2')^-1);

% display the transform
transform

% should just be integer transformations
if (max(max(abs(transform-round(transform))))>0.001)
  disp(sprintf('UHOH: non-integer rotation'));
  return
end
transform = round(transform);

% if this is the identity, then all is good.
if (isequal(transform,eye(3)))
  disp(sprintf('%s and %s are congruent',filename1,filename2));
  return
end

xyswap = 0;
%x = 'x';y = 'y';xdisp = '+x';ydisp = '+y';
x = 1;y = 2;xdisp = '+x';ydisp = '+y';
% should never be in z dimension
if (~isequal(transform(:,3),[0 0 1]'))
  disp(sprintf('UHOH: z dimension transformation?'));
  return
end
% check for x, y swap
if (~isequal(abs(transform(:,2)),[0 1 0]'))
  disp(sprintf('x and y are swapped'));
  xyswap = 1;
  % update transform matrix
  temp = transform(:,2);
  transform(:,2) = transform(:,1);
  transform(:,1) = temp;
end
% check for y flip
if (transform(2,2) == -1)
  disp(sprintf('Y is flipped'));
  y = -2;
  ydisp = '-y';
end
% check for x flip
if (transform(1,1) == -1)
  disp(sprintf('X is flipped'));
  x = -1;
  xdisp = '-x';
end

if (xyswap)
  swapname = sprintf('%s_swapdim%s%s',stripext(filename1),ydisp,xdisp);
%  avwcommand = sprintf('avwswapdim %s %s %s z %s',stripext(filename1),y,x,swapname);
  cbiswapvect = [y x 3];
else
  swapname = sprintf('%s_swapdim%s%s',stripext(filename1),xdisp,ydisp);
%  avwcommand = sprintf('avwswapdim %s %s %s z %s',stripext(filename1),x,y,swapname);
  cbiswapvect = [x y 3];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convert quaternion 2 rotation matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function R = quaternion2rotation(a,b,c,d)

qfac = [1 1 1]';
% assume nifti structure
if (nargin == 1)
  if (isfield(a,'quatern_b'))
    b = a.quatern_b;
    c = a.quatern_c;
    d = a.quatern_d;
    qfac = [1 1 a.pixdim(1)];
    a = real(sqrt(1.0-(b*b+c*c+d*d)));
  % otherwise, passed in length 4 array
  % which has the quaternion values
  else
    b = a(2);c = a(3);d = a(4);a = a(1);
  end
% or, if only three values passed in then
% a is assumed normalized
elseif (nargin == 3)
  d = c;c = b;b = a;
  a = real(sqrt(1.0-(b*b+c*c+d*d)));
end

% formula
R = [ a*a+b*b-c*c-d*d   2*b*c-2*a*d       2*b*d+2*a*c;...
      2*b*c+2*a*d       a*a+c*c-b*b-d*d   2*c*d-2*a*b;...
      2*b*d-2*a*c       2*c*d+2*a*b       a*a+d*d-c*c-b*b ];
%R = (R'*diag(qfac))';
R = (R*diag(qfac));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% getlastdir
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function retval = getlastdir(d)

if (isunix)
  slashes = findstr('/', d);
else
  slashes = findstr('\', d);
end

if (length(slashes) == 0)
  s = 1;
  e = length(d);
elseif (slashes(length(slashes)) == length(d))
  s = slashes(length(slashes)-1)+1;
  e = length(d)-1;
else
  s = slashes(length(slashes))+1;
  e = length(d);
end

retval = d(s:e);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check if this is a file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function retval = isfile(filename)

if (nargin ~= 1)
  help isfile;
  return
end

% open file
fid = fopen(filename,'r');

% check to see if there was an error
if (fid ~= -1)
  fclose(fid);
  retval = 1;
else
  retval = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% stripext
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function retval = stripext(filename)

if (nargin ~= 1)
  help stripext;
  return
end

retval = filename;
dotloc = findstr(filename,'.');
if length(dotloc) > 0
  retval = filename(1:dotloc(length(dotloc))-1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% getext
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function retval = getext(filename)

if (nargin ~= 1)
  help getext;
  return
end

retval = '';
dotloc = findstr(filename,'.');
if (length(dotloc) > 0) && (dotloc(length(dotloc)) ~= length(filename))
  retval = filename(dotloc(length(dotloc))+1:length(filename));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% select an img from a directory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function retval = getimgfromdir(dirname)

retval = [];
% check arguments
if (nargin ~= 1)
  help getimgfromdir;
  return
end

if (isdir(dirname))
  % get files in directory
  thisdir = dir(sprintf('%s/*.img',dirname));
  % if there is only one .img file then this must be it
  if (length(thisdir) == 1)
    retval = sprintf('%s/%s',dirname,thisdir(1).name);
  % if there is more then let the user select 
  elseif (length(thisdir) > 1)
    % display the files
    for i = 1:length(thisdir)
      disp(sprintf('%i: %s',i,thisdir(i).name));
    end
    % get response
    validresponse = 0;
    disp('---------------');
    while (~validresponse)
      r = round(input(sprintf('Select file (%i-%i) to use (0 for none)? ',1,length(thisdir))));
      % response must be valid
      if (length(r) ~= 1) || (r < 0) || (r > length(thisdir))
	continue
      else
	validresponse = 1;
      end
    end
    if (r == 0)
      return
    else
      retval = sprintf('%s/%s',dirname,thisdir(r).name);
    end
  end
else
  disp(sprintf('UHOH: No directory %s',dirname));
  return
end

  
