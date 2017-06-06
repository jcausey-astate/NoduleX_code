"""
Get the (x,y,z) coordinates of each tumor annotated in the supplied 
XML annotation files for the corresponding LIDC-format dicom file.

Usage: 
    {0} DICOM-DIRECTORY [XML-DIRECTORY]

        DICOM-DIRECTORY must contain the dicom image slices (but a single image can be used)
        XML-DIRECTORY   must contain xml files with annotations corresponding to 
                        the dicom images in DICOM-DIRECTORY (or a single xml file can be used)
                        If this parameter is omitted, an XML file must be located
                        in DICOM-DIRECTORY.
        -t N            threshold for malignancy (1-5), omit to see all including "non-nodule" annotations
        -f              file name only (not full path info) in output
        -H              include column headings at top of output
        -m              Matrix-style output (only the numbers: x, y, isNodule, inclusion, malignancy)
"""

from __future__ import print_function
import sys, dicom, re, os
import xml.etree.ElementTree as ET

def usage():
    print(__doc__.format(sys.argv[0]))

def coords_to_list(coords):
    coord_list = []
    for pair in coords:
        coord_list.extend([pair['x'], pair['y']])
    return coord_list

def area_for_polygon(polygon):
    """
        Get the area of a polygon
        Code from: http://stackoverflow.com/a/14115494
    """
    result = 0
    imax = len(polygon) - 1
    for i in range(0,imax):
        result += (polygon[i]['x'] * polygon[i+1]['y']) - (polygon[i+1]['x'] * polygon[i]['y'])
    result += (polygon[imax]['x'] * polygon[0]['y']) - (polygon[0]['x'] * polygon[imax]['y'])
    return abs(result / 2)

def centroid_for_polygon(polygon, non_negative=True):
    """ 
        Get the centroid for a polygon
        Code from: http://stackoverflow.com/a/14115494
    """
    # Special case:  if the annotation only has 1 point, it is already the center:
    if len(polygon) == 1:
        return polygon[0]
    elif len(polygon) < 4:
        # If there are 2 or 3 points, it is really just a line -- return the midpoint.
        return {'x': (polygon[0]['x']+polygon[1]['x']) / 2.0, \
                'y': (polygon[0]['y']+polygon[1]['y']) / 2.0}

    area = area_for_polygon(polygon)
    imax = len(polygon) - 1    
    result_x = 0
    result_y = 0
    if area != 0:
        for i in range(0,imax):
            result_x += (polygon[i]['x'] + polygon[i+1]['x']) * ((polygon[i]['x'] * polygon[i+1]['y']) - (polygon[i+1]['x'] * polygon[i]['y']))
            result_y += (polygon[i]['y'] + polygon[i+1]['y']) * ((polygon[i]['x'] * polygon[i+1]['y']) - (polygon[i+1]['x'] * polygon[i]['y']))
        result_x += (polygon[imax]['x'] + polygon[0]['x']) * ((polygon[imax]['x'] * polygon[0]['y']) - (polygon[0]['x'] * polygon[imax]['y']))
        result_y += (polygon[imax]['y'] + polygon[0]['y']) * ((polygon[imax]['x'] * polygon[0]['y']) - (polygon[0]['x'] * polygon[imax]['y']))
    
        result_x /= (area * 6.0)
        result_y /= (area * 6.0)
    else:
        # area is zero when you have a line or line segments or several coordinates at one pixel.
        # Compute a simple midpoint.
        x_values = [p['x'] for p in polygon]
        y_values = [p['y'] for p in polygon]
        min_x, max_x = min(x_values), max(x_values)
        min_y, max_y = min(y_values), max(y_values)
        result_x = (max_x + min_x) / 2.0
        result_y = (max_y + min_y) / 2.0
    if non_negative:
        result_x = abs(result_x)
        result_y = abs(result_y)
    return {'x': result_x, 'y': result_y}

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def missing_sop_uid(info):
    eprint('Missing SOP UID field in XML file for annotation: {0}'.format(info))
    return "unknown"

# SOPInstanceUID  is same as XML  imageSOP_UID  and corresponds to orig filename

# XML has:
# LidcReadMessage
# ResponseHeader 
#    - SeriesInstanceUid  (maps to directory name, SeriesInstanceUID in DICOM)
# readingSession [...]
#    - annotationVersion
#    - servicingRadiologistID
#    - unblindedReadNodule
#       - noduleID
#       - characteristics
#           - subtlety
#           - internalStructure
#           - calcification
#           - sphericity
#           - margin
#           - lobulation
#           - spiculation
#           - texture
#           - malignancy
#       - roi [...]
#           - imageZposition  (maps to SliceLocation in DICOM)
#           - imageSOP_UID    (maps to SOPInstanceUID in DICOM)
#           - inclusion (bool)
#           - edgeMap [...]
#               - xCoord
#               - yCoord

def get_tumor_centroids(dicomdir, xmldir=None, malignancy_min=None, file_name_only=False, return_file_list=False):
    files       = []
    if os.path.isdir(dicomdir):
        for rt, dr, fl in os.walk(dicomdir):
            files.extend(fl)
            break
    else:
        files    = [dicomdir]
        dicomdir = os.path.abspath(os.path.join(dicomdir, os.pardir))
        if xmldir == None:
            for rt, dr, fl in os.walk(dicomdir):
                for fname in fl:
                    if fname[-4:] == '.xml':
                        if len(files) == 2:
                            eprint("Ambiguous xml files found.")
                        files.append(os.path.join(dicomdir, fname))
                break

    dicomfiles = [f for f in files if f[-4:] == '.dcm']
    xmlfiles   = [f for f in files if f[-4:] == '.xml']

    if xmldir != None:
        if os.path.isdir(xmldir):
            for rt, dr, fl in os.walk(dicomdir):
                if fl[-4:] == '.xml':
                    xmlfiles.extend(fl)
                break    
        else:
            xmlfiles = [xmldir]
    
    if len(xmlfiles) < 1:
        eprint("No XML annotation file found.\n\n")
        usage()
        exit(2)

    nodulesByFile = {}
    namespace     = ""
    dirUID        = ""

    for xfile in xmlfiles:
        tree = ET.parse(os.path.join(dicomdir, xfile))
        root = tree.getroot()
        ns   = ""
        try:
            ns   = re.search('(\{[^}]+)\}', root.tag).group(0)
        except:
            eprint("WARNING: XML seems to be missing namespace -- root tag is {0}".format(root.tag))
        if dirUID == "":
            dirUID = root.find('{0}ResponseHeader/{0}SeriesInstanceUid'.format(ns)).text

        # For each readingSession
        for session in root.iter('{0}readingSession'.format(ns)):
            for nodule in session.iter('{0}unblindedReadNodule'.format(ns)):
                noduleID = nodule.find('{0}noduleID'.format(ns))
                if noduleID != None:
                    noduleID = noduleID.text
                else:
                    noduleID = "unknown"
                characteristics = nodule.find('{0}characteristics'.format(ns))
                malignancy      = None
                if characteristics != None:
                    malignancy      = characteristics.find('{0}malignancy'.format(ns))
                malignancy          = int(malignancy.text) if malignancy != None else None
                if (malignancy_min == None) or ((malignancy != None) and (int(malignancy_min) <= int(malignancy))):
                    for roi in nodule.iter('{0}roi'.format(ns)):
                        coords = []
                        zCoord = roi.find('{0}imageZposition'.format(ns))
                        zCoord = zCoord.text if zCoord != None else 0
                        zCoord = str(float(zCoord))
                        inclusion = roi.find('{0}inclusion'.format(ns))
                        inclusion = inclusion.text if inclusion != None else "FALSE"
                        inclusion = 1 if inclusion.lower() == "true" else 0
                        sopUID = roi.find('{0}imageSOP_UID'.format(ns))
                        sopUID = sopUID.text if sopUID != None else missing_sop_uid("z={0}, nodule={1}".format(zCoord, nID))
                        if not sopUID in nodulesByFile:
                            nodulesByFile[sopUID] = []
                        for edge in roi.iter('{0}edgeMap'.format(ns)):
                            (xCoord, yCoord) = (edge.find('{0}xCoord'.format(ns)).text, edge.find('{0}yCoord'.format(ns)).text)
                            coords.append({'x':float(xCoord), 'y':float(yCoord)})
                        center = centroid_for_polygon(coords)
                        if area_for_polygon(coords) == 0 and len(coords) > 3:
                            eprint("WARNING: Annotation with zero area found; nodule ID {0}, z={1}, # of annotation points: {2}".format(noduleID, zCoord, len(coords)))
                        center['z']          = float(zCoord)
                        center['sopUID']     = sopUID
                        center['noduleID']   = noduleID
                        center['isNodule']   = 1
                        center['inclusion']  = inclusion
                        center['malignancy'] = malignancy if malignancy != None else 0
                        nodulesByFile[sopUID].append(center)
            if malignancy_min == None:
                for nonNodule in session.iter('{0}nonNodule'.format(ns)):
                    nonNoduleID = nonNodule.find('{0}nonNoduleID'.format(ns))
                    if nonNoduleID != None:
                        nonNoduleID = nonNoduleID.text
                    else:
                        nonNoduleID = "_Non-nodule unknown"
                    center = {}
                    zCoord = nonNodule.find('{0}imageZposition'.format(ns))
                    zCoord = zCoord.text if zCoord != None else 0
                    zCoord = str(float(zCoord))
                    sopUID = nonNodule.find('{0}imageSOP_UID'.format(ns))
                    sopUID = sopUID.text if sopUID != None else missing_sop_uid("z={0}, nodule={1}".format(zCoord, nID))
                    if not sopUID in nodulesByFile:
                        nodulesByFile[sopUID] = []
                    locus  = nonNodule.find('{0}locus'.format(ns))
                    if locus == None:
                        continue    
                    (xCoord, yCoord) = (locus.find('{0}xCoord'.format(ns)).text, locus.find('{0}yCoord'.format(ns)).text)
                    center['x'] = float(xCoord)
                    center['y'] = float(yCoord)
                    center['z'] = float(zCoord)
                    center['sopUID']     = sopUID
                    center['noduleID']   = nonNoduleID
                    center['isNodule']   = 0
                    center['inclusion']  = 0
                    center['malignancy'] = 0
                    nodulesByFile[sopUID].append(center)

    result_info = []
    file_list   = []
    # Now, for every dicom file, we can match metadata and output info.
    for dfile in dicomfiles:
        ddat      = dicom.read_file(os.path.join(dicomdir, dfile))
        zCoord    = str(ddat.SliceLocation)
        sopUID    = ddat.SOPInstanceUID
        seriesUID = ddat.SeriesInstanceUID
        sliceNo   = ddat.InstanceNumber
        file_list.append({'sliceNo':sliceNo, 'fileName':os.path.join(dicomdir,dfile), 'seriesUID':seriesUID})
        if sopUID in nodulesByFile:
            for center in nodulesByFile[sopUID]:
                file_name_info = os.path.join(dicomdir,dfile) if not file_name_only else dfile
                center['fileName']  = file_name_info
                center['sliceNo']   = sliceNo
                center['seriesUID'] = seriesUID
                result_info.append(center)
    return result_info if not return_file_list else (result_info, file_list)

def get_dicomdir_files(dicom_dir, include_xml = False, xmldir=None):
    files       = []
    dicom_dir = os.path.normpath(dicom_dir)
    if os.path.isdir(dicom_dir):
        for rt, dr, fl in os.walk(dicom_dir):
            files.extend(fl)
            break
    else:
        files     = [dicom_dir]
        dicom_dir = os.path.abspath(os.path.join(dicom_dir, os.pardir))
        if xmldir == None:
            for rt, dr, fl in os.walk(dicom_dir):
                for fname in fl:
                    if fname[-4:] == '.xml':
                        if len(files) == 2:
                            eprint("Ambiguous xml files found.")
                        files.append(os.path.join(dicom_dir, fname))
                break

    dicomfiles = [f for f in files if f[-4:] == '.dcm']
    xmlfiles   = [f for f in files if f[-4:] == '.xml']
    return (dicomfiles, xmlfiles) if include_xml else dicomfiles


def get_tumor_annotations(dicomdir, xmldir=None, file_name_only=False, return_file_list=False):
    dicomfiles, xmlfiles = get_dicomdir_files(dicomdir, include_xml=True, xmldir=xmldir)
    
    if len(xmlfiles) < 1:
        eprint("No XML annotation file was found.\n\n")
        usage()
        exit(2)

    nodulesByFile = {}
    namespace     = ""
    dirUID        = ""

    for xfile in xmlfiles:
        # print("Reading xml: {}".format(os.path.join(dicomdir, xfile)))
        tree = ET.parse(os.path.join(dicomdir, xfile))
        root = tree.getroot()
        ns   = ""
        try:
            ns   = re.search('(\{[^}]+)\}', root.tag).group(0)
        except:
            eprint("WARNING: XML seems to be missing namespace -- root tag is {0}".format(root.tag))
        if dirUID == "":
            dirUID = root.find('{0}ResponseHeader/{0}SeriesInstanceUid'.format(ns)).text

        # For each readingSession
        for session in root.iter('{0}readingSession'.format(ns)):
            # print("session .. ")
            for nodule in session.iter('{0}unblindedReadNodule'.format(ns)):
                noduleID = nodule.find('{0}noduleID'.format(ns))
                if noduleID != None:
                    noduleID = noduleID.text
                else:
                    noduleID = "unknown"
                characteristics = nodule.find('{0}characteristics'.format(ns))
                malignancy      = None
                if characteristics != None:
                    malignancy      = characteristics.find('{0}malignancy'.format(ns))
                malignancy          = int(malignancy.text) if malignancy != None else None
                for roi in nodule.iter('{0}roi'.format(ns)):
                    # print("nodule .. ")
                    coords = []
                    zCoord = roi.find('{0}imageZposition'.format(ns))
                    zCoord = zCoord.text if zCoord != None else 0
                    zCoord = str(float(zCoord))
                    inclusion = roi.find('{0}inclusion'.format(ns))
                    inclusion = inclusion.text if inclusion != None else "FALSE"
                    inclusion = 1 if inclusion.lower() == "true" else 0
                    sopUID = roi.find('{0}imageSOP_UID'.format(ns))
                    sopUID = sopUID.text if sopUID != None else missing_sop_uid("z={0}, nodule={1}".format(zCoord, nID))
                    if not sopUID in nodulesByFile:
                        nodulesByFile[sopUID] = []
                    for edge in roi.iter('{0}edgeMap'.format(ns)):
                        (xCoord, yCoord) = (edge.find('{0}xCoord'.format(ns)).text, edge.find('{0}yCoord'.format(ns)).text)
                        coords.append({'x':float(xCoord), 'y':float(yCoord)})
                    center = centroid_for_polygon(coords)
                    center['boundary'] = coords
                    center['coords']   = coords_to_list(coords)
                    if area_for_polygon(coords) == 0 and len(coords) > 3:
                        eprint("WARNING: Annotation with zero area found; nodule ID {0}, z={1}, # of annotation points: {2}".format(noduleID, zCoord, len(coords)))
                    center['z']          = float(zCoord)
                    center['sopUID']     = sopUID
                    center['noduleID']   = noduleID
                    center['isNodule']   = 1
                    center['inclusion']  = inclusion
                    center['malignancy'] = malignancy if malignancy != None else 0
                    # print('center: {}'.format(center))
                    nodulesByFile[sopUID].append(center)
            
            for nonNodule in session.iter('{0}nonNodule'.format(ns)):
                # print("non nodule .. ")
                nonNoduleID = nonNodule.find('{0}nonNoduleID'.format(ns))
                if nonNoduleID != None:
                    nonNoduleID = nonNoduleID.text
                else:
                    nonNoduleID = "_Non-nodule unknown"
                center = {}
                zCoord = nonNodule.find('{0}imageZposition'.format(ns))
                zCoord = zCoord.text if zCoord != None else 0
                zCoord = str(float(zCoord))
                sopUID = nonNodule.find('{0}imageSOP_UID'.format(ns))
                sopUID = sopUID.text if sopUID != None else missing_sop_uid("z={0}, nodule={1}".format(zCoord, nID))
                if not sopUID in nodulesByFile:
                    nodulesByFile[sopUID] = []
                locus  = nonNodule.find('{0}locus'.format(ns))
                if locus == None:
                    continue    
                (xCoord, yCoord) = (locus.find('{0}xCoord'.format(ns)).text, locus.find('{0}yCoord'.format(ns)).text)
                center['x'] = float(xCoord)
                center['y'] = float(yCoord)
                center['z'] = float(zCoord)
                center['sopUID']     = sopUID
                center['noduleID']   = nonNoduleID
                center['isNodule']   = 0
                center['inclusion']  = 0
                center['malignancy'] = 0
                nodulesByFile[sopUID].append(center)

    result_info = []
    file_list   = []
    # Now, for every dicom file, we can match metadata and output info.
    # print('for dfile in {}'.format(dicomfiles))
    for dfile in dicomfiles:
        # print('reading {}'.format(os.path.join(dicomdir, dfile)))
        ddat      = dicom.read_file(os.path.join(dicomdir, dfile))
        zCoord    = str(ddat.SliceLocation)
        sopUID    = ddat.SOPInstanceUID
        seriesUID = ddat.SeriesInstanceUID
        sliceNo   = ddat.InstanceNumber
        patient   = os.path.basename(dicomdir)
        file_list.append({'sliceNo':sliceNo, 'fileName':os.path.join(dicomdir,dfile), 'seriesUID':seriesUID, 'patient':patient, 'dicom_directory':dicomdir})
        if sopUID in nodulesByFile:
            # print("matched sopUID {}".format(sopUID))
            for center in nodulesByFile[sopUID]:
                file_name_info = os.path.join(dicomdir,dfile) if not file_name_only else dfile
                center['fileName']  = file_name_info
                center['filePath']  = os.path.join(dicomdir, dfile)
                center['sliceNo']   = sliceNo
                center['seriesUID'] = seriesUID
                center['patient']   = patient
                center['dicom_directory']  = dicomdir
                # print("appending center at slice {}".format(sliceNo))
                result_info.append(center)
        else:
            # print("sopUID didn't match any from xml: {}".format(sopUID))
            pass
    return result_info if not return_file_list else (result_info, file_list)

def print_results(res, dicomdir=None, numeric_only=False):
    if len(res) == 0:
        location = "{0}: ".format(dicomdir) if dicomdir != None else ""
        eprint("{0}No tumor annotations found corresponding to input dataset.".format(location))
    else:
        for info in res:
            if not numeric_only:
                print('\t'.join([str(x) for x in [info['fileName'], info['noduleID'], info['sliceNo'], info['x'], info['y'], info['z'], info['isNodule'], info['inclusion'], info['malignancy']]]))
            else:
                print('\t'.join([str(x) for x in [info['sliceNo'], info['x'], info['y'], info['isNodule'], info['inclusion'], info['malignancy']]]))

def run_standalone():
    if len(sys.argv) < 2:
        usage()
        sys.exit(1)
    mal_thresh      = None
    xmlfile         = ""
    dicomdir        = sys.argv[1]
    xmldir          = None if len(sys.argv) < 3 else sys.argv[2]
    file_name_only  = False
    file_header     = False
    matrix_style    = False
    token = xmldir
    if token != None and token[0] == '-':
        xmldir = None
    i = 2 
    while i < len(sys.argv):
        token = sys.argv[i]
        if token == '-t':
            mal_thresh = int(sys.argv[i+1])
            i += 1
        elif token == '-f':
            file_name_only = True
        elif token == '-H':
            file_header = True
        elif token == '-m':
            matrix_style = True
        i += 1
    res = get_tumor_centroids(dicomdir, xmldir, malignancy_min=mal_thresh, file_name_only=file_name_only)
    if file_header:
        if matrix_style:
            print('\t'.join(["sliceNo", "x", "y", "z", "isNodule", "inclusion", "malignancy"]))
        else:    
            print('\t'.join(["fileName", "noduleID", "sliceNo", "x", "y", "z", "isNodule", "inclusion", "malignancy"]))
    print_results(res, dicomdir, matrix_style)

if __name__ == "__main__":
    run_standalone()
