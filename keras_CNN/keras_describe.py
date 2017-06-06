"""
Describe a keras model (load it and create graph and description output)

keras_describe.py  model_file

"""


from __future__ import print_function
from keras.preprocessing.image import ImageDataGenerator
from keras.models import Sequential
from keras.layers import Dense, Dropout, Activation, Flatten
from keras.layers import Convolution2D, MaxPooling2D
from keras.optimizers import SGD, Adagrad, Adadelta, RMSprop
from keras.utils import np_utils
try:
    from keras.utils.visualize_util import plot
except ImportError:
    from keras.utils.vis_utils import plot_model as plot
from keras.models import load_model
from keras.models import model_from_json
from keras.callbacks import EarlyStopping, TensorBoard
import sys, random, pickle, numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
from datetime import date


filename = sys.argv[1]

model = None
if filename[-3:3] in ['hd5', 'hdf', 'df5']:
	model = load_model(filename)
else:
	with open(filename, 'rU') as json_file:
		model = model_from_json(json_file.read())

print('Network Layout:')
model.summary()

plot_name = '.'.join(filename.split('.')[0:-1]) + ".png"

plot(model, to_file=plot_name, show_shapes=True)

