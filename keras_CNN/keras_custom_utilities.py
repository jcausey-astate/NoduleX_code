"""
Provides some functions for convenience that are not directly available in 
Keras framework.
"""

def crop(img, shape):
    """
    Crop img to shape (which can potentially be either smaller than img,
    or larger - in which case img will be zero-padded).
    @param img      Array of 4 dimensions (images, z, y, x) where images are
                    all images of the original shape (z, y, x).
    @param shape    Array of 3 dimensions (z, y, x) representing the goal
                    (cropped) shape for each image.
    """
    sz, sy, sx = shape
    hy, hx     = int(sy / 2), int(sx / 2)
    # print("Cropping from {0} to {1}".format(img.shape, shape))
    n_chan, mz, my, mx = img.shape
    cy, cx     = int(my / 2), int(mx / 2)
    if sz == mz and sy == my and sx == mx:
        # print("No crop needed.")
        return np.copy(img)
    results    = []
    for i in range(n_chan):
        cropped = np.zeros(shape).astype(img[i].dtype)
        cx      = int(max(0, sx-mx) / 2)
        cy      = int(max(0, sy-my) / 2)
        cz      = int(max(0, sz-mz) / 2)
        ox      = int(max(0, mx-sx) / 2)
        oy      = int(max(0, my-sy) / 2)
        oz      = int(max(0, mz-sz) / 2)
        # print("cz {0}, cy {1}, cx {2}, oz {3}, oy {4}, ox {5}".format(cz, cy, cx, oz, oy, ox))
        cropped[cz:sz-cz, cy:sy-cy, cx:sx-cx] = img[i][oz:mz-oz, oy:my-oy, ox:my-ox].copy()
        results.append(cropped)
    results = np.array(results)
    # print("Finished cropping. Returning results shape {0}".format(results.shape))
    return results