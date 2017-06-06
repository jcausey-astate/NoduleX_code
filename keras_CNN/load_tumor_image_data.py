"""
Load lung-ct image data for training a CNN, where the data was created in picklefile 
or hdf5 file formed by the 'create_data_files_from_tumors.py' script
"""

import numpy as np, random, h5py, os

def file_ext(filename):
    return os.path.splitext(filename)[1][1:] if filename is not None else None

def load_all_data(filename=None, normalize=False, malignancy_to_class=None, window_normalize=False):
    if file_ext(filename).lower() in ['hdf5', 'hd5', 'hdf', 'h5']:
        return load_all_from_hdf5(filename, normalize=normalize, malignancy_to_class=malignancy_to_class, window_normalize=window_normalize)
    else:
        return load_all_from_picklefile(filename)

def load_train_test_data(filename=None, test_pct=0.25, neg_bias=None, shuffle_seed=None, normalize=False, malignancy_to_class=None, window_normalize=False):
    neg_bias = 0.5 if neg_bias is None else neg_bias
    if file_ext(filename) in ['hdf5', 'hd5', 'hdf', 'h5']:
        return load_train_test_from_hdf5(filename, test_pct, neg_bias, shuffle_seed=shuffle_seed, normalize=normalize, malignancy_to_class=malignancy_to_class, window_normalize=window_normalize)
    else:
        if normalize or window_normalize:
            raise RuntimeError("normalize not supported from legacy picklefiles...")
        return load_train_test_from_picklefile(filename, test_pct, neg_bias, shuffle_seed=shuffle_seed, malignancy_to_class=malignancy_to_class)

def load_metadata(filename=None):
    if not file_ext(filename).lower() in ['hdf5', 'hd5', 'hdf', 'h5']:
        raise RuntimeError("This operation is only supported for HDF5 files.")
    import h5py, numpy as np
    fp = h5py.File(filename, 'r')
    hd5_metadata = {}
    avoid_list   = ['nodule_images', 'nodule_classes']
    for key in fp.keys():
        if key not in avoid_list:
            hd5_metadata[key] = fp[key].value
    return hd5_metadata        

def load_all_from_picklefile(filename=None, test_pct=0.25, neg_bias=None, shuffle_seed=None, malignancy_to_class=None):
    """
    @param filename str     name of the file to load
    @return (X,y)           all data from the file
    """
    neg_bias = 0.5 if neg_bias is None else neg_bias
    import pickle, numpy as np
    if filename == None:
        filename = 'LIDC-IDRI.pickle'
    print("Loading data from {0}".format(filename))
    X,y = pickle.load(open(filename, 'rb'))
    return (np.array(X), np.array(y))

def load_train_test_from_picklefile(filename=None, test_pct=0.25, neg_bias=None, shuffle_seed=None, malignancy_to_class=None):
    """
    @param filename str     name of the file to load
    @param test_pct float   the percentage of examples to use for testing
    @param neg_bias float   [optional] the percentage of overall examples that should be negative
                            (or None to use the data as-is)  Default = None
    @return ((Train_X, Train_Y), (Test_X, Test_Y)) train and test set tuples are returned
    """
    import pickle, numpy as np, random
    if filename == None:
        filename = 'LIDC-IDRI.pickle'
    neg_bias = 0.5 if neg_bias is None else neg_bias
    print("Loading data from {0}".format(filename))
    X,y = pickle.load(open(filename, 'rb'))
    n_examples = len(X)
    if shuffle_seed != None:
        random.seed(shuffle_seed)

    if neg_bias != None:
        negatives = [i for i in range(n_examples) if y[i] == [0]]
        neg_count = len(negatives)
        pos_count = n_examples - neg_count
        neg_goal  = int(round(neg_bias * n_examples))
        pos_goal  = n_examples - neg_goal
        # print("Before: neg count: {0}; goal: {1} - pos count: {2}; goal: {3}".format(neg_count, neg_goal, pos_count, pos_goal))
        
        # Remove either some postive examples or some negative examples until we get the correct
        # negative-to-positive bias (if requested)
        victims = []
        removed = 0
        if neg_count > neg_goal: 
            random.shuffle(negatives)
            while neg_count > neg_goal:
                victims.append(negatives.pop())
                removed   += 1
                neg_count -= 1
                neg_goal   = int(round(neg_bias * (n_examples - removed)))
        elif pos_count > pos_goal:
            positives = [i for i in range(n_examples) if y[i] != 0]
            random.shuffle(positives)
            while pos_count > pos_goal:
                victims.append(positives.pop())
                removed   += 1
                pos_count -= 1
                neg_goal   = int(round(neg_bias * (n_examples - removed)))
                pos_goal   = n_examples - removed - neg_goal    
        X_new = []
        y_new = []
        for i in range(n_examples):
            if not i in victims:
                X_new.append(X[i])
                y_new.append(y[i])
        X = X_new
        y = y_new
        n_examples = len(X)
    pos_goal = neg_count + pos_count - neg_goal
    # print("After: neg count: {0}; goal: {1} - pos count: {2}; goal: {3}".format(neg_count, neg_goal, pos_count, pos_goal))

    test_examples = int(round(test_pct * n_examples))
    indices = range(0,n_examples)
    random.shuffle(indices)
    X_test  = [X[i] for i in indices[0:test_examples]]
    y_test  = [y[i] for i in indices[0:test_examples]]
    X_train = [X[i] for i in indices[test_examples:]]
    y_train = [y[i] for i in indices[test_examples:]]
        
    print("Divided input into {0} training and {1} test examples.".format(len(X_train), len(X_test)))
    return (np.array(X_train), np.array(y_train)), (np.array(X_test), np.array(y_test))

def load_all_from_hdf5(filename=None,  normalize=False, malignancy_to_class=None, window_normalize=False):
    """
    @param filename str     name of the file to load
    @return (X,y)           all data from the file
    """
    import h5py, numpy as np
    if filename == None:
        filename = 'LIDC-IDRI.hd5'
    print("Loading data from {0}".format(filename))
    fp = h5py.File(filename, 'r')
    y  = fp['nodule_classes'].value
    if malignancy_to_class is not None:
        if len(malignancy_to_class) != 6:
            raise Exception("malignancy_class mapping must contain exactly 6 values, one for each malignancy level 0 - 5")
        mal = fp['nodule_malignancy']          
        for i in range(len(y)):
            y[i] = [malignancy_to_class[int(mal[i])]]
        y = [[v] for v in y]
    X = fp['nodule_images'].value
    X = np.array(X, dtype='float32')

    # Normalize if requested
    if normalize or window_normalize:
        Xmin = fp['nodule_pixel_min'].value
        Xmax = fp['nodule_pixel_max'].value
        print("{}Normalizing...".format("Window " if window_normalize else ""))
        for idx, img in enumerate(X):
            nXmin = Xmin[idx] if not window_normalize else -1000.0
            nXmax = Xmax[idx] if not window_normalize else  4096.0
            print("               {} {}".format(nXmin, nXmax))
            X[idx] = (np.array(X[idx]) - nXmin) / (nXmax - nXmin)
            X[idx][X[idx] < 0] = 0
            X[idx][X[idx] > 1] = 1

    # Remove any cases where class is < 0:
    y = np.array(y)
    X = X[np.where(y.ravel() >= 0)]
    y = y[np.where(y.ravel() >= 0)]
    y = np.array(y)
    return (X, y)

def load_train_test_from_hdf5(filename=None, test_pct=0.25, neg_bias=None, batch_size=64, shuffle_seed=None, normalize=False, malignancy_to_class=None, window_normalize=False):
    """
    @param filename str     name of the file to load
    @param test_pct float   the percentage of examples to use for testing
    @param neg_bias float   [optional] the percentage of overall examples that should be negative
                            (or None to use the data as-is)  Default = None
    @return ((Train_X, Train_Y), (Test_X, Test_Y)) train and test set tuples are returned
    """
    import h5py, numpy as np, random
    if filename == None:
        filename = 'LIDC-IDRI.hd5'
    neg_bias = 0.5 if neg_bias is None else neg_bias
    print("Loading data from {0}".format(filename))
    if shuffle_seed != None:
        random.seed(shuffle_seed)
    # Use the y (class labels) vector to get length and class count stats; loading the 
    # X values (the image cubes) could be memory-prohibitive.  Instead, use an iterable
    # generator object to stream batches that meet the class distribution requirements
    # on demand.
    tts = TestTrainSet(filename, test_pct, neg_bias, batch_size, normalize=normalize, malignancy_to_class=malignancy_to_class, window_normalize=window_normalize)
    return tts.get_train_test()


class AutoCrop(object):
    def __init__(self, images, shape):
        self._needs_crop     = True
        self._len            = len(images)
        self._orig_images    = images
        self._images         = np.nditer(images)
        self._out_shape      = shape
        sz, sy, sx           = shape
        self.shape           = (images.shape[0], sz, sy, sx)
        self._orig_shape     = images.shape
        n_frames, mz, my, mx = self._orig_shape
        self._n_frames       = n_frames
        self._dtype          = images[0].dtype
        if sz == mz and sy == my and sx == mx:
            self._needs_crop = False
        cx            = int(max(0, sx-mx) / 2)
        cy            = int(max(0, sy-my) / 2)
        cz            = int(max(0, sz-mz) / 2)
        ox            = int(max(0, mx-sx) / 2)
        oy            = int(max(0, my-sy) / 2)
        oz            = int(max(0, mz-sz) / 2)
        self._center  = (cz, cy, cx)
        self._ocenter = (oz, oy, ox)
        
    def __iter__(self):
        return self

    def _crop_many(self, imgs):
        results  = []
        for slice in imgs:
            results.append(self._crop_one(slice))
        return np.array(results, self._dtype)

    def _crop_one(self, img):
        result = None
        if self._needs_crop:
            sz, sy, sx      = self._out_shape
            mp, mz, my, mx  = self._orig_shape
            cz, cy, cx      = self._center
            oz, oy, ox      = self._ocenter        
            cropped         = np.zeros(self._out_shape).astype(self._dtype)
            # print("cz {0}, cy {1}, cx {2}, oz {3}, oy {4}, ox {5}".format(cz, cy, cx, oz, oy, ox))
            cropped[cz:sz-cz, cy:sy-cy, cx:sx-cx] = img[oz:mz-oz, oy:my-oy, ox:my-ox]
            result = cropped
        else:
            result = img
        return result

    def _crop(self, img):
        result = None
        if len(img.shape) > 3:
            result = self._crop_many(img)
        else:
            result = self._crop_one(img)
        return result

    def __next__(self):
        return self._crop(next(self._images))
    
    def next(self):
        return self.__next__()

    def reset(self):
        self._images.reset()

    def __len__(self):
        return self._len

    def __getitem__(self, index):
        return self._crop(self._orig_images[index])

    def __copy__(self):
        return AutoCrop(self._orig_images.copy(), self._out_shape)

    def __deepcopy__(self):
        return self.__copy__()



class TestTrainSet(object):
    def __init__(self, hdf_filename, test_pct=0.25, neg_bias=0.5, batch_size=64,  normalize=False, malignancy_to_class=None, window_normalize=False):       
        neg_bias = 0.5 if neg_bias is None else neg_bias
        self._hdf_filename        = hdf_filename
        self._neg_bias            = neg_bias
        self._test_pct            = test_pct
        self._batch_size          = batch_size
        self._test_location       = 0
        self._train_location      = 0
        self._test_indices        = []
        self._train_indices       = []
        self._malignancy_to_class = malignancy_to_class
        self._normalize           = normalize
        self._Xmin                = None
        self._Xmax                = None
        self._window_normalize    = window_normalize
        if malignancy_to_class is not None and len(malignancy_to_class) != 6:
            raise Exception("malignancy_class mapping must contain exactly 6 values, one for each malignancy level 0 - 5")
        # Open the hdf file 
        self._hdf_file = h5py.File(self._hdf_filename, 'r')
        # Get info on classes and makeup of the dataset by examining the y values (classes):
        y          = self._hdf_file['nodule_classes'].value        
        if self._malignancy_to_class is not None:
            if malignancy_to_class is not None:
                mal = self._hdf_file['nodule_malignancy']          
                for i in range(len(y)):
                    y[i] = [malignancy_to_class[int(mal[i])]]
        if self._normalize:
            self._Xmin = self._hdf_file['nodule_pixel_min']
            self._Xmax = self._hdf_file['nodule_pixel_max']    
        n_examples = len(y)
        negatives  = [i for i in range(n_examples) if y[i] == [0]]
        positives  = [i for i in range(n_examples) if y[i][0] > 0]
        neg_count  = len(negatives)
        pos_count  = len(positives)
        n_examples = neg_count + pos_count
        neg_goal   = int(min(neg_count, round(neg_bias * n_examples)))
        pos_goal   = n_examples - neg_goal
        if pos_goal > pos_count:
            neg_goal = int(round(pos_count * (1-neg_bias+0.5)))
            pos_goal = pos_count
        # print("Before: neg count: {0}; goal: {1} - pos count: {2}; goal: {3}".format(neg_count, neg_goal, pos_count, pos_goal))
        # randomly choose neg_goal negatives and pos_goal positives:
        selected_indices = list(np.random.choice(negatives, size=(min(neg_count,neg_goal)), replace=False))
        selected_indices.extend(list(np.random.choice(positives, size=(min(pos_count, pos_goal)), replace=False)))
        # print("n examples: {0}".format(len(selected_indices)))
        n_examples = len(selected_indices)
        np.random.shuffle(selected_indices)
        test_examples = int(round(test_pct * n_examples))
        # print("test_examples: {0}".format(test_examples))
        self._test_indices  = selected_indices[0:test_examples]
        self._train_indices = selected_indices[test_examples:]

    def get_train_test_batch(self):
        train_indices = self._train_indices[self._train_location:self._train_location + self._batch_size]
        test_indices  = self._test_indices[self._test_location:self._test_location + self._batch_size]
        self._train_location += self._batch_size
        self._test_location  += self._batch_size
        h5f     = self._hdf_file
        h5f_X   = h5f['nodule_images']
        h5f_y   = h5f['nodule_classes'] if self._malignancy_to_class is None else h5f['nodule_malignancy']
        X_train = np.take(h5f_X, train_indices, axis=0)
        y_train = np.take(h5f_y, train_indices, axis=0)
        X_test  = np.take(h5f_X, test_indices,  axis=0)
        y_test  = np.take(h5f_y, test_indices,  axis=0)
        Xmin_train, Xmax_trian = None, None
        Xmin_test,  Xmax_test  = None, None
        if self._normalize or self._window_normalize:
            Xmin_train = np.take(h5f['nodule_pixel_min'], train_indices, axis=0)
            Xmax_train = np.take(h5f['nodule_pixel_max'], train_indices, axis=0)
            Xmin_test  = np.take(h5f['nodule_pixel_min'], test_indices,  axis=0)
            Xmax_test  = np.take(h5f['nodule_pixel_max'], test_indices,  axis=0)
            print("{}Normalizing...".format("Window " if self._window_normalize else ""))
            for idx, img in X_train:
                # print("               {} {}".format(Xmin_train[idx], Xmax_train[idx]))
                nXmin = Xmin_train[idx] if not self._window_normalize else -1000.0
                nXmax = Xmax_train[idx] if not self._window_normalize else  4096.0
                X_train[idx] = (np.array(X_train[idx]) - nXmin) / (nXmax - nXmin)
                X_train[idx][X_train[idx] < 0] = 0
                X_train[idx][X_train[idx] > 1] = 1 
            for idx, img in X_test:
                # print("               {} {}".format(Xmin_test[idx], Xmax_test[idx]))
                nXmin = Xmin_test[idx] if not self._window_normalize else -1000.0
                nXmax = Xmax_test[idx] if not self._window_normalize else  4096.0
                X_test[idx]  = (np.array(X_test[idx])  - nXmin)  / (nXmax - nXmin)
                X_test[idx][X_test[idx] < 0] = 0
                X_test[idx][X_test[idx] > 1] = 1 
            # print("min-max now: tr {} {} ts {} {}".format(X_train.min(), X_train.max(), X_test.min(), X_test.max()))

        if self._malignancy_to_class is not None:
            y_train = [[int(self._malignancy_to_class[v[0]])] for v in y_train]
            y_test  = [[int(self._malignancy_to_class[v[0]])] for v in y_test]
        
        # print("y_train: {}".format(y_train))
        y_train = np.array(y_train)
        y_test  = np.array(y_test)

        return (X_train, y_train), (X_test, y_test)

    def get_train_test(self):
        train_indices = self._train_indices
        test_indices  = self._test_indices
        self._train_location = 0
        self._test_location  = 0
        h5f     = self._hdf_file
        h5f_X   = h5f['nodule_images']
        h5f_y   = h5f['nodule_classes'] if self._malignancy_to_class is None else h5f['nodule_malignancy']
        X_train = np.take(h5f_X, train_indices, axis=0)
        y_train = np.take(h5f_y, train_indices, axis=0)
        X_test  = np.take(h5f_X, test_indices,  axis=0)
        y_test  = np.take(h5f_y, test_indices,  axis=0)
        if self._normalize or self._window_normalize:
            Xmin_train = np.take(h5f['nodule_pixel_min'], train_indices, axis=0)
            Xmax_train = np.take(h5f['nodule_pixel_max'], train_indices, axis=0)
            Xmin_test  = np.take(h5f['nodule_pixel_min'], test_indices,  axis=0)
            Xmax_test  = np.take(h5f['nodule_pixel_max'], test_indices,  axis=0)
            print("{}Normalizing...".format("Window " if self._window_normalize else ""))
            for idx, img in enumerate(X_train):
                # print("               {} {}".format(Xmin_train[idx], Xmax_train[idx]))
                nXmin = Xmin_train[idx] if not self._window_normalize else -1000.0
                nXmax = Xmax_train[idx] if not self._window_normalize else  4096.0
                X_train[idx] = (np.array(X_train[idx]) - nXmin) / (nXmax - nXmin)
                X_train[idx][X_train[idx] < 0] = 0
                X_train[idx][X_train[idx] > 1] = 1 
            for idx, img in enumerate(X_test):
                # print("               {} {}".format(Xmin_test[idx], Xmax_test[idx]))
                nXmin = Xmin_test[idx] if not self._window_normalize else -1000.0
                nXmax = Xmax_test[idx] if not self._window_normalize else  4096.0
                X_test[idx]  = (np.array(X_test[idx])  - nXmin)  / (nXmax - nXmin)
                X_test[idx][X_test[idx] < 0] = 0
                X_test[idx][X_test[idx] > 1] = 1 
            # print("min-max now: tr {} {} ts {} {}".format(X_train.min(), X_train.max(), X_test.min(), X_test.max()))

        if self._malignancy_to_class is not None:
            y_train = [[int(self._malignancy_to_class[v[0]])] for v in y_train]
            y_test  = [[int(self._malignancy_to_class[v[0]])] for v in y_test]

        # print("y_train: {}".format(y_train))
        y_train = np.array(y_train)
        y_test  = np.array(y_test)

        return (X_train, y_train), (X_test, y_test)

    def __iter__(self):
        return self

    def __next__(self):
        if self._train_location >= len(self._train_indices)  or  self._test_location >= len(self._test_indices):
            # Reset:
            np.shuffle(self._train_indices)
            np.shuffle(self._test_indices)
            self._train_location = 0
            self._test_location  = 0
            raise StopIteration()
        return self.get_train_test_batch()
    next = __next__

def normalize_data(X):
    return normalize_individual(X)

def normalize_set_wise(X):
    X  = X.astype('float32')
    X /= X.max()
    return X

def normalize_individual(X):
    X = X.astype(np.float32)
    for idx, x in enumerate(X):
        x -= x.min()
        x /= x.max()
        X[idx] = x
    return X
