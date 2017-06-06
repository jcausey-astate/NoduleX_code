"""
Evaluate a Keras model on a fresh dataset.
"""

from __future__ import print_function
from keras.models import Model
from keras.optimizers import SGD, Adagrad, Adadelta, RMSprop
from keras.utils import np_utils
from keras.models import load_model
from keras.models import model_from_json
from keras import backend as K
import sys, os, numpy as np
from load_tumor_image_data import *

nb_classes    = 2
batch_size    = 64

def main():
    args = get_args()
    model_file   = args['model_file']
    data_file    = args['data_file']
    weights_file = None
    split_model  = False
    if args['model_weights'] != '':
        weights_file = args['model_weights']
        split_model = True
    elif os.path.splitext(model_file)[1] == '.json':
        weights_file = os.path.splitext(model_file)[0] + '.weights.hd5'
        split_model = True
    window_normalize = args['window_normalize']
   
    # Give some feedback on settings
    if not args['normalize'] and (window_normalize == False):
        print("Using raw images.")
    if window_normalize:
        print("Using window normalization.")
    # Load the data file
    (X,y) = load_all_data(data_file, normalize=args['normalize'], window_normalize=window_normalize)
    # Feedback on the data
    print("X shape: {} ; X[0] shape: {}  X[0][2] shape: {}".format(X.shape, X[0].shape, X[0][2].shape))
    img_rows, img_cols = X[0][2].shape
    print('img_rows: {0}, img_cols: {1}'.format(img_rows, img_cols))
    print('X shape:', X.shape)
    print('{0} samples ({1} (+), {2} (-))'.format(X.shape[0], y.sum(), len(y)-y.sum()))

    # convert class vectors to binary class matrices
    ground_truth = y
    y = np_utils.to_categorical(y, nb_classes)

    # Check validity of model file:
    model = None
    if split_model:
        print('Loading model design {0}\nand weights {1}'.format(model_file, weights_file))
        with open(model_file, 'rU') as json_file:
            model = model_from_json(json_file.read())
        model.load_weights(weights_file)
    else:
        print('Loading full model file {0}'.format(model_file))
        model = load_model(model_file)
    
    if args['show_layout']:
        print('Network Layout:')
        model.summary()
        print('\n\n')
    
    y_pred = model.predict(X, batch_size=batch_size, verbose=1)
    result = quantify_results(y_pred[:,1], ground_truth, auc=True)

    if args['list_cases']:
        metadata = load_metadata(data_file)
        print('')
        print_case_table(y_pred, ground_truth, metadata)
        print('\n\n')

    if args['expression_file'] != None:
        print('\nSaving expression vectors in {0}.'.format(args['expression_file']))
        import csv
        feature_layer = find_feature_layer(model, index=int(args['expression_layer_index']))
        print("Getting output at layer index {0}. Layer output shape: {1}".format(feature_layer, model.layers[feature_layer].output_shape))
        expression    = get_output_from_layer(model, X, layer_index=feature_layer)
        with open(args['expression_file'], 'wb') as csvfile:
            writer = csv.writer(csvfile, delimiter='\t', quoting=csv.QUOTE_MINIMAL)
            for vec in expression:
                writer.writerow(vec)

    if args['predictions_file'] != None:
        # predictions are just expression at the last layer, reduced to a single floating-point value if
        # the number of classes is 2.
        print('\nSaving prediction values in {0}.'.format(args['predictions_file']))
        import csv
        with open(args['predictions_file'], 'wb') as csvfile:
            writer = csv.writer(csvfile, delimiter='\t', quoting=csv.QUOTE_MINIMAL)
            for vec in y_pred:
                try:
                    if len(vec) == 2:
                        vec = [vec[1]]
                except:
                    vec = [vec]
                writer.writerow(vec)        

    print('(+) - Total: {0} ; True: {1} ; False {2}'.format(result['npos'],result['tp'],result['fp']))
    print('(-) - Total: {0} ; True: {1} ; False {2}'.format(result['nneg'],result['tn'],result['fn']))
    print('auc: {0} acc: {1}, sensitivity: {2}, specificity: {3}, precision: {4}'.format(result['auc'], result['acc'], result['sens'],result['spec'], result['prec']))
    return 0

def get_output_from_layer(model, X, layer_index=None, train_mode=False, n_categories=2):
    """
    get the output from a specific layer given the input; defaults to the last
    layer before a reduction to classes (<=n_categories)
    """
    mode = 0 if not train_mode else 1
    if layer_index == None:
        layer_index = find_feature_layer(model, n_categories=n_categories)
     
    get_nth_layer_output = K.function([model.layers[0].input, K.learning_phase()],
                                      [model.layers[layer_index].output])

    layer_output = get_nth_layer_output([X, mode])[0]

    return layer_output


def find_last_feature_layer(model, n_categories=2):
    unwanted_layers  = ['Dropout', 'Pooling']
    last_layer_index = len(model.layers) - 1
    last_layer       = model.layers[last_layer_index]
    while last_layer.output_shape[-1] <= n_categories \
          or any(ltype in str(type(last_layer)) for ltype in unwanted_layers):
          last_layer_index -= 1
          last_layer = model.layers[last_layer_index]
    return last_layer_index

def find_feature_layer(model, index=-2):
    unwanted_layers     = ['Dropout', 'Pooling']
    direction           = np.sign(index)
    i_orig              = index
    index               = abs(index)
    current_layer_index = len(model.layers) - 1
    current_layer       = model.layers[current_layer_index]
    if direction >= 0:
        current_layer       = model.layers[0]
        current_layer_index = 0
    while (direction != 0 and index >= 0) or any(ltype in str(type(current_layer)) for ltype in unwanted_layers):
          current_layer_index += direction
          current_layer = model.layers[current_layer_index]
          index -= 1
    return current_layer_index


def quantify_results(predictions, ground_truth, auc=False):
    fp = fn = tp = tn = npos = nneg = 0
    total = len(predictions)
    if auc:
        from sklearn import metrics
        auc_value   = round(metrics.roc_auc_score(np.ravel(ground_truth), np.ravel(predictions)), 3)
        predictions = np.round(predictions)

    for idx, predicted in enumerate(predictions):
        if is_positive(ground_truth[idx]):
            npos += 1
        else:
            nneg += 1
        if is_positive(predicted) == is_positive(ground_truth[idx]):
            if is_positive(ground_truth[idx]):
                tp += 1
            else:
                tn += 1
        else:
            if is_positive(predicted):
                fp += 1
            else:
                fn += 1
    epsilon = 1e-20 # used to avoid divide-by-zero in rare cases where all example are seen as a single class
    acc  = round((tp + tn) / float(total), 3)
    sens = round(tp  / (float(npos)    + epsilon), 3)
    spec = round(tn  / (float(nneg)    + epsilon), 3)
    ppv  = round(tp  / (float(tp + fp) + epsilon), 3)
    npv  = round(tn  / (float(tn + fn) + epsilon), 3)
    result = {
        'fp': fp,
        'tp': tp,
        'fn': fn,
        'tn': tn,
        'total': total,
        'tot': total,
        'npos': npos,
        'nneg': nneg,
        'acc': acc,
        'sens': sens,
        'tpr': sens,
        'spec': spec,
        'spc': spec,
        'tnr': spec,
        'ppv': ppv,
        'prec': ppv,
        'npv': npv
    }
    if auc:
        result['auc'] = auc_value
    return result

def print_case_table(y_pred, ground_truth, metadata={}):
    meta_keys = metadata.keys()
    headings  = ['true class', 'predicted', 'correct?']
    headings.extend(meta_keys)
    print('\t'.join(headings))
    for idx, prediction in enumerate(y_pred):
        results = [str(x) for x in [ground_truth[idx][0], prediction, str(ground_truth[idx][0]==prediction)]]
        for key in meta_keys:
            try:
                results.append(str(metadata[key][idx][0]))
            except TypeError:
                results.append(str(metadata[key][idx]))
        print('\t'.join(results))

def array_like(x):
    import collections
    return isinstance(x, collections.Sequence) or isinstance(x, np.ndarray)

def is_positive(cls):
    if array_like(cls) and len(cls) > 1:
        return cls[0] == 0
    else:
        return cls != 0 if not array_like(cls) else cls[0] != 0

def get_args():
    import argparse
    # construct the argument parser and parse the arguments
    ap = argparse.ArgumentParser(prog='{0}'.format(os.path.basename(sys.argv[0])))
    ap.add_argument("--list", dest="list_cases", action='store_true', help="List results case-by-case (verbose output).")
    ap.add_argument("--layout", dest="show_layout", action='store_true', help="Show network layout.")
    ap.add_argument("model_file", metavar='model-file' , help = "Keras model file (.hd5 or .json)")
    ap.add_argument("model_weights", metavar='model-weights', nargs='?', default="", 
        help = "Keras weights file (.hd5); optional if model file was .hd5 or if name is same as model file except for extension.")
    ap.add_argument("data_file", metavar='data-file' , help = "HDF5 file containing dataset to evaluate.")
    ap.add_argument("--raw", dest='normalize', action='store_false', help="Use raw images; no normalization.")
    ap.add_argument("--window", dest='window_normalize', action='store_true', help="Perform HU window normalization.")
    ap.add_argument("-x", dest = 'expression_file', required = False, 
        help = "If given, \"expression levels\" for each item in data file are saved here in CSV-compatible format.")
    ap.add_argument("-l", dest = 'expression_layer_index', required = False, default = -2,
        help = """Used in conjunction with '-x' to select a specific layer's expression output; positive values 
                  start from the beginning, negative values from the end.  Default is last layer before reducing 
                  to classes (-2); pooling and dropout layers do not count.""")
    ap.add_argument("-p", dest = 'predictions_file', required=False, 
        help = "If given, \"predictions\" (floating-point) for each item in data file are saved here in CSV-compatible format.")
    args = vars(ap.parse_args())
    return args

if __name__ == '__main__':
    # construct the argument parser and parse the arguments
    sys.exit(main())