# NoduleX
Supporting code for the paper _"Highly accurate model for prediction of lung nodule malignancy with CT scans"_.

## Instructions
After cloning or downloading this repository, extract the files from (http://bioinformatics.astate.edu/NoduleX/NoduleX_data.tar.gz) into the `data` directory.  

Many of the scripts contained here have several command-line options available.  Run the scripts with the `--help` option to see a usage listing.

### Requirements
* Python2.7, `pip`
    - A requirements file `NoduleX_python_requirements.txt` is provided listing required Python packages.  You can install them using:
        - `pip install -r NoduleX_python_requirements.txt`
    - Setting up a virtual environment is recommended.
* QIF feature extraction requires Octave (version tested was 4.2.0), or MATLAB (with some modifications to the helper scripts, see QIF_extraction/README.md).
* A POSIX-compatible system (Linux, Mac OS, or Linux shell under Windows) is assumed; many scripts given are written in Bash shell syntax.

### Running the CNN models against validation data
Use the script `keras_CNN/keras_evaluate.py`, providing the correct model .json file (from `data/CNN_models`) and matching weights .hd5 file (from `data/CNN_weights`) and a dataset .hd5 file (from `data/CNN_datasets`).  For example:

```bash
python keras_CNN/keras_evaluate.py \
    --window \
    data/CNN_models/CNN47.json \
    data/CNN_weights/12v45_weights.hd5 \
    data/CNN_datasets/S12vS45_s47_VALIDATION.hd5
```

### Training CNN models with training data
Use the script `keras_CNN/keras_retrain_model.py`, providing the correct model .json file (from `data/CNN_models`) and a training dataset .hd5 file (from `data/CNN_datasets`).  For example:

```bash
python keras_CNN/keras_retrain_model.py \
    --window \
    -s 47 \
    data/CNN_models/CNN47.json \
    data/CNN_datasets/S12vS45_s47_TRAIN.hd5 \
    /tmp/CNN47_retrain_checkpoint_dir
```

The final model weights will be saved in the working directory (as a .hd5 file), and checkpoints will be placed in the directory `/tmp/CNN47_retrain_checkpoint_dir` (you can customize this of course).

Depending on the number of epochs and batch size you choose (default is 200 and 64), the model may overfit.  Examine the checkpoint models as well as the final model to determine the best overall performance.  (Typically the best will be one of the last 3 checkpoints or the final model.)  Training is stochastic, so repeated training will yield different results.

### Building datasets from LIDC-IDRI
Start by downloading the data files for LIDC-IDRI (
https://wiki.cancerimagingarchive.net/display/Public/LIDC-IDRI) and extract the DOI folder in `data`.

Run the script `dicom_and_image_tools/simplify_doi_structure.sh` against the DOI direcotry (you may choose to create symlinks with the `-s` option).  This produces a directory structure that is flattened with naming based on patient identifiers.  It could also be useful to run `dicom_and_image_tools/rename_dicom_by_position.py` against each of the patient directories to name the .dcm files themselves in order of 'sliceNo' (not necessary, but it makes the files easier to reason about).

#### Extracting CNN Cubes
Extract input volumes for the CNN by using the script `dicom_and_image_tools/create_data_file_from_seed_points.py`.  Candidate nodule lists are in the directory `data/nodule_lists`.  Run the script using the top-level of your flattened data directory from the previous step as `dicom_dir`.

#### Creating Segmentation Masks for QIF Feature Computation
The QIF feature extraction code requires the original image ("grey image") and an image where all pixels representing the ROI are set to 1 while all others are set to 0 ("binary image").  To create the binary images, use `dicom_and_image_tools/segment_to_binary_image.py` as follows:

* **For "nodules"**:  For any nodule with a malignancy rating 1-5, segmentations are provided by LIDC-IDRI.  Run the `segment_to_binary_image.py` script with the `--candidates` option pointing to the candidates file (from `data/nodule_lists`) and the `--segmented-only` option.  This works for the "S12vS45", "S1vS45" and "NvNN_nodule-only" candidate lists.
* **For "non-nodules"**: For the "non-nodule" dataset (the "NvNN_non-nodule-only" candidate list), a segmentation must be algorithmically generated for each "non-nodule" seed point.  Run the `segment_to_binary_image.py` script with the `--candidates` option pointing to the candidates file (from `data/nodule_lists`).  The segmentation process will take some time.

#### Converting to Analyze format for QIF Feature Computation
The QIF Feature extraction code requires Analyze format for its input files; LIDC-IDRI data is in DICOM format.  To convert (both the "grey" and "binary", see above) DICOM files to Analyze, the tool `dicom_and_image_tools/dicom_to_analyze.py` is provided.  Run it for each patient's scan (producing the "grey" images), and for each nodule segmentation ("binary image") file.

For example, if you placed your binary DICOM images in a directory `binary_dicom`, you could something similar to the following to convert all nodules for a single patient:

```bash
p=<YOUR-LIDC-IDRI-PATIENT-ID->; \
for n in `ls -d data/binary_dicom/$p/*` ; do \
    echo "Converting nodule $n" ; \
    python dicom_and_image_tools/dicom_to_analyze.py \
        "$n" \
        "data/binary_analyze/$p/$(basename $n)" \
        && echo "OK" \
        || echo "FAILED converting nodule $n" \
;done
```

Where `<YOUR-LIDC-IDRI-PATIENT-ID>` is the patient ID for the patient whose nodules you are converting.  Adjust paths according to your local directory layout, as necessary.

#### Extracting QIF Features
See the file `README.md` in the `QIF_extraction` directory for the steps required to compute QIF features given the "grey" and "binary" images in Analyze format as described above.





