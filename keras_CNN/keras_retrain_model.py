'''
Re-train a keras model.

Run with --help for usage info.
'''

from __future__ import print_function
from keras.preprocessing.image import ImageDataGenerator
from keras.models import Sequential
from keras.layers import Dense, Dropout, Activation, Flatten
from keras.layers import Convolution2D, MaxPooling2D, ZeroPadding2D
from keras.optimizers import SGD, Adagrad, Adadelta, RMSprop
from keras.utils import np_utils
from keras.models import load_model, model_from_json
from keras.regularizers import l2, activity_l2, l1, activity_l1, l1l2, activity_l1l2
from keras.callbacks import EarlyStopping, TensorBoard, ModelCheckpoint
import sys, os, random, pickle, numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
import datetime
from load_tumor_image_data import *
from image_5ch import *

batch_size        = 64
nb_classes        = 2
nb_epoch          = 200
data_augmentation = True
test_fraction     = 0.2
neg_bias          = 0.5

# we use 5 slices for each "cube"
img_channels = 5

def load_data(filename=None, test_pct=0.25, neg_bias=None, normalize=True, window_normalize=False):
    return load_train_test_data(filename, test_pct, neg_bias, normalize=normalize, window_normalize=window_normalize)

def main():
    args            = get_args()
    img_rows        = args['size']
    img_cols        = img_rows
    img_shape       = (img_channels, img_rows, img_cols)
    nb_epoch        = args['epochs']
    batch_size      = args['batch_size']
    is_cropped      = False
    model_file      = args['model_file']
    orig_model_file = model_file
    data_file       = args['data_file']
    chkpt_prefix    = '/tmp/keras_checkpoint'
    if args['chkpt_prefix'] != None:
        chkpt_prefix = args['chkpt_prefix']
    chkpt_prefix = chkpt_prefix[0:-1] if chkpt_prefix[-1] == '/' else chkpt_prefix
    if not os.path.isdir(chkpt_prefix):
        os.makedirs(chkpt_prefix)
    
    model_file = ""
    chkpt_file = ""
    if(data_file != None):
        model_file += '.'.join(data_file.split('.')[0:-1]) + "."
    model_basename = os.path.basename(model_file)
    chkpt_prefix   = os.path.join(chkpt_prefix, model_basename)
    datestr        = "{0}".format(datetime.datetime.now().isoformat().split('.')[0:-1][0])
    chkpt_file     = chkpt_prefix + "keras_model_" + datestr + "_ep-{epoch:02d}_acc-{val_acc:.3f}_loss-{val_loss:.3f}.hd5"

    # the data, shuffled and split between train and test sets
    do_normalize = not args['raw']
    if not do_normalize:
        print("Using raw input images.")
   
    (X_train, y_train), (X_test, y_test) = load_data(data_file, test_fraction, neg_bias=neg_bias, normalize=do_normalize, window_normalize=args['window'])
    orig_shape           = X_train[0].shape
    orig_channels, orig_rows, orig_cols = orig_shape
    print('Input data img_rows: {0}, img_cols: {1}'.format(orig_rows, orig_cols))
    print('X_train shape:', X_train.shape)
    if orig_shape != img_shape:
        print('Cropping data to (x,y,z): ({0},{1},{2}) from ({3},{4},{5})'.format(img_cols, img_rows, img_channels, orig_cols, orig_rows, orig_channels))
        is_cropped = True
        chkpt_file = chkpt_prefix + "keras_model_{0}x{1}_".format(img_cols, img_rows) + datestr + "_ep-{epoch:02d}_acc-{val_acc:.3f}_loss-{val_loss:.3f}.hd5"
    print('{0} training samples ({1} (+), {2} (-))'.format(X_train.shape[0], y_train.sum(), len(y_train)-y_train.sum()))
    print('{0} testing samples ({1} (+), {2} (-))'.format(X_test.shape[0], y_test.sum(), len(y_test)-y_test.sum()))
    max_rank = int(max(img_rows, img_cols))

    # convert class vectors to binary class matrices
    y_train = np_utils.to_categorical(y_train, nb_classes)
    y_test  = np_utils.to_categorical(y_test, nb_classes)

    split_model  = False
    if os.path.splitext(orig_model_file)[1] == '.json':
        split_model = True

    model = None
    if split_model:
        print('Loading model design {0}'.format(orig_model_file))
        with open(orig_model_file, 'rU') as json_file:
            model = model_from_json(json_file.read())
    else:
        print('Loading full model file {0}'.format(orig_model_file))
        model = load_model(orig_model_file)
        weights = model.get_weights()
    
    opt = Adadelta()
    model.compile(loss='categorical_crossentropy',
                  optimizer=opt,
                  metrics=['accuracy'])

    X_train = X_train.astype('float32')
    X_test  = X_test.astype('float32')
    X_train /= X_train.max()
    X_test  /= X_test.max()

    early_stopping = EarlyStopping(monitor='val_loss', patience=2)
    print("Ckpt file name: " + chkpt_file)
    checkpointing  = ModelCheckpoint(chkpt_file, monitor='val_loss', verbose=0, save_weights_only=False, save_best_only=True)

    hist = None

    if not data_augmentation:
        print('Not using data augmentation.')
        hist = model.fit(X_train, y_train,
                  batch_size=batch_size,
                  nb_epoch=nb_epoch,
                  validation_data=(X_test, y_test),
                  shuffle=True,
                  callbacks=[checkpointing])
    else:
        print('Using real-time data augmentation.')

        # this will do preprocessing and realtime data augmentation
        datagen = ImageDataGenerator(
            featurewise_center=False,  # set input mean to 0 over the dataset
            samplewise_center=False,  # set each sample mean to 0
            featurewise_std_normalization=False,  # divide inputs by std of the dataset
            samplewise_std_normalization=False,  # divide each input by its std
            zca_whitening=False,  # apply ZCA whitening
            rotation_range=180,  # randomly rotate images in the range (degrees, 0 to 180)
            width_shift_range=0.25,  # randomly shift images horizontally (fraction of total width)
            height_shift_range=0.25,  # randomly shift images vertically (fraction of total height)
            horizontal_flip=False,  # randomly flip images
            vertical_flip=False)  # randomly flip images

        # compute quantities required for featurewise normalization
        # (std, mean, and principal components if ZCA whitening is applied)
        datagen.fit(X_train, img_shape)

        # fit the model on the batches generated by datagen.flow()
        hist = model.fit_generator(datagen.flow(X_train, y_train,
                            batch_size=batch_size),
                            samples_per_epoch=X_train.shape[0],
                            nb_epoch=nb_epoch,
                            validation_data=(X_test, y_test),
                            callbacks=[checkpointing])
    cropped_dims = ""
    if is_cropped:
        cropped_dims = "_{0}x{1}".format(img_cols, img_rows)
    model_info       = "{3}_{0}_acc-{1:.3}_loss-{2:.3}".format(datestr, hist.history['val_acc'][-1], hist.history['val_loss'][-1], cropped_dims)
    model_file_base  = os.path.basename(model_file) + "keras_model_{0}".format(model_info)
    model_file       = model_file_base + '.hd5'
    
    model_json  = model_file_base + '.json'
    model_weight_file = model_file_base + '.weights.hd5'
    with open(model_json, 'w') as fp_json:
        fp_json.write(model.to_json())
    model.save_weights(model_weight_file)

    history_file = "keras_train_stats_{0}.pickle".format(model_info)
    pickle.dump(hist.history, open(history_file, 'wb'))


def get_args():
    import argparse
    # construct the argument parser and parse the arguments
    ap = argparse.ArgumentParser(prog='{0}'.format(os.path.basename(sys.argv[0])))
    ap.add_argument("model_file", metavar='model-file' , help = "Keras model file (.hd5 or .json)")
    ap.add_argument("--raw", action='store_true', help = "Use raw images; do not normalize.")
    ap.add_argument("--window", action='store_true', help = "Perform a 'window normalization' so that the range -1000 to 4096HU is scaled to [0,1]")
    ap.add_argument("--epochs", dest='epochs', default = 200, type=int, help = "Number of epochs to train (default 200).")
    ap.add_argument("--batch-size", dest='batch_size', default=64, type=int, help = "Number of samples per batch.")
    ap.add_argument("-s", "--size", dest='size', default=21, type=int, help = "Size (width/height) of input patches.")
    ap.add_argument("data_file", metavar='data-file' , help = "HDF5 file containing dataset to evaluate.")
    ap.add_argument("chkpt_prefix", metavar='checkpoint-dir', nargs='?', default=None, help = "Directory to place checkpoint files.")
    args = vars(ap.parse_args())
    return args

if __name__ == '__main__':
    # construct the argument parser and parse the arguments
    sys.exit(main())