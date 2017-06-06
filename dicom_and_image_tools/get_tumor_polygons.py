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
        -m              Matrix-style output (only the numbers: isNodule, inclusion, malignancy, z, x0, y0, x1, y1, ... , xN, yN)
"""

from __future__ import print_function
import sys, dicom, re, os
import xml.etree.ElementTree as ET

global standalone ; standalone = False

def usage(msg):
    global standalone
    eprint(msg + '\n\n')
    if standalone:
        print(__doc__.format(sys.argv[0]))
    else:
        raise RuntimeError(msg)

def area_for_polygon(polygon):
    """
        Get the area of a polygon
        Code from: http://stackoverflow.com/a/14115494
    """
    result = 0
    imax   = len(polygon) - 1
    for i in range(0,imax):
        result += (polygon[i]['x'] * polygon[i+1]['y']) - (polygon[i+1]['x'] * polygon[i]['y'])
    result += (polygon[imax]['x'] * polygon[0]['y']) - (polygon[0]['x'] * polygon[imax]['y'])
    return abs(result / 2.0)

def centroid_for_polygon(coords, non_negative=True):
    """ 
        Get the centroid for a polygon
        Code from: http://stackoverflow.com/a/14115494
    """
    polygon = []
    for i in range(0, len(coords)-1, 2):
        polygon.append({'x': coords[i], 'y': coords[i+1]})
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

def coords_to_list(coords):
    coord_list = []
    for pair in coords:
        coord_list.extend([pair['x'], pair['y']])
    return coord_list

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

def get_tumor_polygons(dicomdir, xmldir=None, malignancy_min=0, file_name_only=False):
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
        usage("No XML annotation file found.")
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
                malignancy          = int(malignancy.text) if malignancy != None else 0
                if int(malignancy_min) <= int(malignancy):
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
                        area   = area_for_polygon(coords)
                        center = centroid_for_polygon(coords_to_list(coords))
                        if area == 0 and len(coords) > 3:
                            eprint("NOTICE: Annotation with zero area found; nodule ID {0}, z={1}, # of annotation points: {2}".format(noduleID, zCoord, len(coords)))
                        if area > 0:
                            polygon               = {}
                            polygon['coords']     = coords_to_list(coords)
                            polygon['x']          = center['x']
                            polygon['y']          = center['y']
                            polygon['z']          = float(zCoord)
                            polygon['area']       = area
                            polygon['sopUID']     = sopUID
                            polygon['noduleID']   = noduleID
                            polygon['isNodule']   = 1
                            polygon['inclusion']  = inclusion
                            polygon['malignancy'] = malignancy if malignancy != None else 0
                            nodulesByFile[sopUID].append(polygon)

    result_info = [];
    # Now, for every dicom file, we can match metadata and output info.
    for dfile in dicomfiles:
        ddat      = dicom.read_file(os.path.join(dicomdir, dfile))
        zCoord    = str(ddat.SliceLocation)
        sopUID    = ddat.SOPInstanceUID
        seriesUID = ddat.SeriesInstanceUID
        sliceNo   = ddat.InstanceNumber
        if sopUID in nodulesByFile:
            for polygon in nodulesByFile[sopUID]:
                file_name_info = os.path.join(dicomdir,dfile) if not file_name_only else dfile
                polygon['fileName']  = file_name_info
                polygon['sliceNo']   = sliceNo
                polygon['seriesUID'] = seriesUID
                result_info.append(polygon)
    return result_info

def print_results(res, dicomdir=None, numeric_only=False):
    if len(res) == 0:
        location = "{0}: ".format(dicomdir) if dicomdir != None else ""
        eprint("{0}No tumor annotations found corresponding to input dataset.".format(location))
    else:
        for info in res:
            if not numeric_only:
                print('\t'.join([str(x) for x in [info['fileName'], info['noduleID'], info['sliceNo'], info['inclusion'], info['malignancy'], info['z']] + [str(x) for x in info['coords']] ]))
            else:
                print('\t'.join([str(x) for x in [info['sliceNo'], info['inclusion'], info['malignancy']] + [str(x) for x in info['coords']] ]))

def run_standalone():
    if len(sys.argv) < 2:
        usage("Not enough arguments.")
        sys.exit(1)
    mal_thresh      = 0
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
    res = get_tumor_polygons(dicomdir, xmldir, malignancy_min=mal_thresh, file_name_only=file_name_only)
    if file_header:
        if matrix_style:
            print('\t'.join(["sliceNo", "inclusion", "malignancy", "x0", "y0", "...", "...", "xN", "yN"]))
        else:
            print('\t'.join(["fileName", "noduleID", "sliceNo", "inclusion", "malignancy", "z", "x0", "y0", "...", "...", "xN", "yN"]))
    print_results(res, dicomdir, matrix_style)

if __name__ == "__main__":
    standalone = True
    run_standalone()
