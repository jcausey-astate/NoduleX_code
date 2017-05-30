# Extracting QIF features

The script `lna_parallel_n.sh` is provided to simplify extracting features from a batch of nodules.  The requirements and details for using this script are presented below.

## Software requirements:
The script expects Octave (tested with version 4.2.0) with the _statistics_ and _image_ packages installed.  The code itself will also run under MATLAB, but you will need to modify `lna_batch_driver.sh` and `lna_by_file.sh` to change the Octave commands to the corresponding MATLAB commands.

## Directory and file naming conventions
In order for the `lna_parallel_n.sh` script to associate the files correctly, it is important that the following directory structure and naming conventions be used.

* You will need a directory for the original CT scan images in Analyze format (referred to as "grey" images here), and a second directory for segmentation mask images (referred to as "binary" images here).
    - For this example, assume the name `grey_dir` is used for the "grey" images, and `mask_dir` is used for the "binary" segmentation masks.
* Beneath the "grey" directory, create one directory for each patient.
    - For example, in `grey_dir` create a directory for each patient, using a unique patient ID for each, like `grey_dir/PATIENT-0001`, `grey_dir/PATIENT-0002`, etc.
* Beneath the "binary" directory, create one directory for each nodule, where the directory name should contain the patieht ID followed by a tilde (`~`) symbol followed by a nodule ID.
    - For example, in `mask_dir`, if `PATIENT-0001` has two nodules `Nodule1` and `Nodule2`, you would create the directories `mask_dir/PATIENT-0001~Nodule1` and `mask_dir/PATIENT-0001~Nodule2`.

To summarize, patients must all have a unique ID and that their files will be stored in a directory named with this ID.  Nodules must also have a unique ID, and it must be combined with the patient ID in the format `PATIENT-ID~NODULE-ID` (separated by a tilde character); nodule files should be placed into a directory with the combined patient/nodule ID as its name.

## Prepare the "grey" and "binary" images
The extraction code requires that a matching "grey" (original CT scan) and "binary" (segmentation mask) image exists for each nodule to be examined.

The "binary" segmentation mask file must be the same size (in terms of X,Y,Z dimension) as the "grey" file, and must have all voxel values set to zero except for voxels belonging to the nodule, which should all be set to 1.

You can use the same "grey" file for all nodules for a specific patient, but you will need a separate "binary" segmentation mask file for each nodule.  The input files must be in Analyze (.img) format. 

For example, if you have a patient PATIENT-0001 with two nodules Nodule1 and Nodule2, you would need to prepare the following:

* `grey_dir/PATIENT-0001/PATIENT-0001.img , grey_dir/PATIENT-0001/PATIENT-0001.hdr` Analyze format "grey" (original CT image) for PATIENT-0001
* `mask_dir/PATIENT-0001~Nodule1/PATIENT-0001~Nodule1.img , mask_dir/PATIENT-0001~Nodule1/PATIENT-0001~Nodule1.hdr` - Analyze format "binary" segmentation mask for PATIENT-0001, Nodule1.
* `mask_dir/PATIENT-0001~Nodule2/PATIENT-0001~Nodule2.img , mask_dir/PATIENT-0001~Nodule2/PATIENT-0001~Nodule2.hdr` - Analyze format "binary" segmentation mask for PATIENT-0001, Nodule2.

## Create one or more nodule list files
You can run the extraction process in parallel -- the number of parallel instances is determined by the number of nodule list files you prepare (each will be run in parallel with the others).

For example, given the files detailed above for two nodules belonging to PATIENT-0001, you can create two nodule lists, and run two processes in parallel.  Create the following files:

* `nodule-list.1`
* `nodule-list.2`

With the following contents:

**`nodule-list.1`**:

```
PATIENT-0001~Nodule1
```

**`nodule-list.2`**
```
PATIENT-0001~Nodule2
```

## Running in parallel
With the list files prepared as detailed above, you can run as follows:

```
./lna_parallel_n.sh \
    -g grey_dir \
    -b mask_dir \
    -o output_dir \
    nodule-list.1 nodule-list.2
```

Where `grey_dir` and `mask_dir` match the directories where you have stored the "grey" and "binary" files, and `output_dir` is the location where you want the output to be written.

## Convert output to CSV for analysis
The output (.hd5) file is not very user-friendly for downstream analysis.  Scripts are provided to help convert it to a more usable form:

* `batch_matlab_features_to_csv.sh` - Used to extract the features from the HDF5 (.hd5) files into easy-to-use CSV (comma-separated) format.
    - This script will produce one CSV file per nodule, matching the one HDF5 file per nodule produced by the Octave/MATLAB code.

**Example**
```
./batch_matlab_features_to_csv.sh \
    --output_dir /csv_output_dir \
    /output_dir/*.hd5
```

* `collect_features_to_matrix.py` - A Python script used to combine the individual nodule-level CSV files into a single CSV file with one column per nodule.
    - For the "Highly accurate model for prediction of lung nodule malignancy with CT scans" paper, we used the following options:
        + `--remove-rows 26 43 44`

**Example**
```
python collect_features_to_matrix.py    \
    /combined_output_dir/combined2d.csv  \
    /csv_output_dir/*features2d.csv \
    --remove-rows 26 43 44
```

# Author Info
The code in the `Matlab_Source` subdirectory was mainly written by David Politte and David Gierada.  The driver scripts in this directory were mostly written by Jason L Causey; some code in this directory was written by Justin Porter, and the utility _cell2csv.m_ was written by Sylvain Fiedler.  Other code authors and attribution information is listed in header comments in code.
