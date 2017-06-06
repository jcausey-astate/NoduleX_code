"""
Library used by other scripts for internally managing a data structure
of "candidates" (ROIs of either the "centroid" type or the "polygon" type
where polygons have a segmentation associated with them and centroids do
not).  Also provides functions for managing groups of candidates.
"""
import sys, os, numpy as np

def calculate_group_xyz(group, include_volume=False):
    '''
    Given a group of candidates, calculate the common x,y,z location
    of the group.  May also calculate and return the (voxel) volume if 
    `include_volume` is True.
    '''
    cx, cy, cz  = (0,0,0)
    volume      = 0
    members     = group
    if 'members' in group:
        members = group['members']
    for slc in members:
        area    = float(slc['area']) if float(slc['area']) > 0 else 1
        cx     += float(slc['x']) * area
        cy     += float(slc['y']) * area
        cz     += float(slc['sliceNo']) * area
        volume += area
    if volume == 0:
        volume = float(len(members)) if len(members) > 0 else 1
    cx /= volume
    cy /= volume
    cz /= volume
    return (cx,cy,cz) if not include_volume else (cx,cy,cz,volume)

def get_group_nodule_flags(group):
    '''
    Gets the three nodule-level metadata items (in LIDC-IDRI, these 
    correspond to "is-nodule", "inclusion", and "malignancy", but they may
    be used arbitrarily for other datasets).
    '''
    is_nodule   = 0
    inclusion   = 0
    malignancy  = 0
    members     = group
    if 'members' in group:
        members = group['members']
    for slc in members:
        malignancy += int(slc['malignancy'])
        is_nodule   = is_nodule | int(slc['isNodule'])
        inclusion   = inclusion | int(slc['inclusion'])
    malignancy = int(round(malignancy / float(len(members))))  # Malignancy is the average of all opinions (rounded to nearest int)
    return is_nodule, inclusion, malignancy

def get_candidate_seriesUID(group):
    '''
    Retrieve the 'seriesUID' attribute for a single candidate or a 
    group of candidates.
    '''
    seriesUID = group['seriesUID'] if 'seriesUID' in group else None
    if seriesUID is None and 'members' in group:
        group = group['members']
    if seriesUID is None and type(group) == type([]):
        for candidate in group:
            seriesUID = get_candidate_seriesUID(candidate)
            if seriesUID is not None:
                break # All it takes is one...
    return seriesUID

def set_group_stats(group):
    '''
    Calculate the identifying metadata for a group of candidates and 
    generate an ID for the group.
    '''
    cx, cy, cz, volume = calculate_group_xyz(group, include_volume=True)
    is_nodule, inclusion, malignancy = get_group_nodule_flags(group)
    group['isNodule']   = is_nodule
    group['inclusion']  = inclusion
    group['malignancy'] = malignancy
    group['x'] = cx
    group['y'] = cy
    group['z'] = cz
    group['sliceNo']  = int(round(cz))
    group['volume']   = volume
    group['noduleID'] = make_unique_id_string(cx, cy, cz, is_nodule, inclusion, malignancy)
    return group

def make_unique_id_string(x, y, z, is_nodule, inclusion, malignancy, prefix=None):
    '''
    Given the information required for the unique ID, put it into the common 
    string format and return the ID string.
    '''
    prefix = "{}~".format(prefix) if prefix is not None else ""
    if z is not None:
        n_id = "{}{}-{}-{}~{}-{}-{}".format(prefix, int(round(x)), int(round(y)), int(round(z)), is_nodule, inclusion, malignancy)
    else:
        n_id = "{}{}-{}~{}-{}-{}".format(prefix, int(round(float(x))), int(round(float(y))), is_nodule, inclusion, malignancy)
    return n_id
        
def nodule_unique_id(candidate, prefix=None):
    '''
    Create a unique ID string for a candidate based on its X,Y,slceNo location 
    and either a provided `prefix` or a hash of the 'seriesUID' attribute.
    '''
    n_id = ""
    if 'members' in candidate:
        candidate = candidate['members']
    if type(candidate) == type([]):
        cx, cy, cz  = calculate_group_xyz(candidate)
        is_nodule, inclusion, malignancy = get_group_nodule_flags(candidate)
        seriesUID   = get_candidate_seriesUID(candidate)
        if seriesUID is not None and prefix is None:
            prefix = hash_string(seriesUID)
        n_id = make_unique_id_string(cx, cy, cz, is_nodule, inclusion, malignancy, prefix)
    else:
        x          = candidate['x']
        y          = candidate['y']
        is_nodule  = candidate['isNodule']   if 'isNodule'   in candidate else 0
        inclusion  = candidate['inclusion']  if 'inclusion'  in candidate else 0
        malignancy = candidate['malignancy'] if 'malignancy' in candidate else 0
        seriesUID  = candidate['seriesUID']  if 'seriesUID'  in candidate else None
        if seriesUID is not None and prefix is None:
            prefix = hash_string(seriesUID)
        z = None
        if 'sliceNo' in candidate:
            z = candidate['sliceNo']
        n_id = make_unique_id_string(x, y, z, is_nodule, inclusion, malignancy, prefix)
    return n_id

def img_hash(img):
    '''
    Compute the sha1 hash of the contents of an (np array) image.
    '''
    import hashlib

    sha = hashlib.sha1()
    sha.update(img.data)
    return sha.hexdigest()

def hash_string(value):
    '''
    Compute the sha1 hash of the contents of an (np array) image.
    '''
    import hashlib

    sha = hashlib.sha1()
    sha.update(value)
    return sha.hexdigest()

def fill_file_paths(candidates, dicom_dir):
    '''
    Populate the 'filePath' attribute of all candidates, relative to `dicom_dir`.
    '''
    import dicom, dicom_dir_utilities
    single_candidate = False
    if not type(candidates) == type([]):
        single_candidate = True
        candidates = [candidates]
    
    file_list = dicom_dir_utilities.get_dicomdir_files(dicom_dir, full_path=True)
    files_by_slice = {}
    
    for fname in file_list:
        slice_no = int(dicom.read_file(fname).InstanceNumber)
        files_by_slice[slice_no] = fname

    for idx, candidate in enumerate(candidates):
        c_slice = int(candidate['sliceNo'])
        candidates[idx]['filePath'] = files_by_slice[c_slice]

    return candidates if not single_candidate else candidates[0]

if __name__ == '__main__':
    print("This file is designed to be included as a library only.")
    sys.exit(1)