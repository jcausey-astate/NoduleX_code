% CBI - NIFTI tools for Matlab
%
%     (c) Copyright Jonas Larsson 2005-2008
% 
%     This file is part of cbiNifti: an Nifti-1 I/O library for Matlab and Octave.
% 
%     cbiNifti is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     cbiNifti is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with cbiNifti.  If not, see <http://www.gnu.org/licenses/>.
%
%
% NOTE: these tools require the mlpipe library for compressed file support.
%
% TOOLS FOR READING NIFTI FILES
% =============================
% cbiReadNiftiHeader.m
% - Reads in the header of a Nifti-1 file
% cbiReadNifti.m
% - Reads in the data and header of a Nifti-1 file
% 
% TOOLS FOR WRITING NIFTI FILES
% =============================
% cbiWriteNiftiHeader.m
% - Writes a Nifti-1 header
% cbiWriteNifti.m
% - Writes a Nifti-1 file (data and header)
%
% TOOLS FOR CONVERTING NIFTI FILES
% ================================
% cbiSwapNiftiDimensions.m
% - Swaps the dimensions of a Nifti data set or file. Analogous to FSL avwswapdim, but correctly sets qform & sform data.
%
% NIFTI UTILITIES
% ===============
% cbiCreateNiftiHeader.m 
% - Creates a new header and/or checks an existing header for consistency and valid fields
% cbiHomogeneousToQuaternion.m 
% - Converts a 4x4 or 3x3 matrix to quaternions
% cbiQuaternionToHomogeneous.m
% - Converts quaternions and qoffset data into a 4x4 homogeneous matrix
% cbiSetNiftiQform.m
% - Sets the qform and quaternion information of a Nifti header
% cbiSetNiftiSform.m
% - Sets the sform and srow information of a Nifti header
% cbiMatlabDatatype2Nifti.m
% - Converts matlab data type in string format to a Nifti-1 integer code; see nifti1.h for details 
% cbiNiftiDatatype2Matlab.m
% - Converts a Nifti-1 integer data tyoe code into matlab data type in string format
% cbiSizeofNifti.m
% - Returns the size in bytes per data point of a given Nifti-1 format.
% cbiParseNiftiFilename.m
% - Parses file name into image, header, extension, and compression fragments.
% nifti1.h
% - Header file of the nifti-1 library written by Robert Cox; describes the Nifti-1 format
