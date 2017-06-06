import unittest
import sitkstrats
import SimpleITK as sitk  # pylint: disable=F0401
import numpy as np

# pylint: disable=missing-docstring
# pylint: disable=invalid-name

class TestConsensusSeg(unittest.TestCase):
    '''test sitkstrats.segmentation_union with various fudged datasets'''

    def test_orthogonal_segs(self):
        a1 = np.zeros((20, 20, 20), dtype='uint8')
        a1[3:8, 3:8, 3:8] = 1

        a2 = np.zeros(a1.shape, dtype='uint8')
        a2[13:17, 13:17, 13:17] = 1

        (i1, i2) = [sitk.GetImageFromArray(a) for a in (a1, a2)]

        # segmentation_union([i1, i2], )


class TestCOMCalc(unittest.TestCase):
    '''test sitkstrats.com_calc with various fudged data'''

    def setUp(self):
        self.MAX_SIZE = 1e10  # pylint: disable=C0103
        self.MIN_SIZE = 0  # pylint: disable=C0103

        self.test = np.zeros((10, 10, 10), dtype='int64')
        self.test[2:4, 2:4, 2:4] = 1

        self.ones_img = sitk.GetImageFromArray(
            np.ones(self.test.shape, dtype='int64'))


    def test_simple(self):
        self.test[2:4, 2:4, 2:4] = 1
        img = sitk.GetImageFromArray(self.test)

        (seeds, info) = sitkstrats.com_calc(
            img, self.MAX_SIZE, self.MIN_SIZE, self.ones_img)

        self.assertEqual(len(seeds), 1)
        self.assertEqual(len(seeds[0]), 3)
        self.assertEqual(type(seeds[0][0]), int)
        self.assertEqual(seeds[0], [2, 2, 2])
        self.assertEqual(info['nseeds'], 1)
        self.assertEqual(info['min_size'], self.MIN_SIZE)
        self.assertEqual(info['max_size'], self.MAX_SIZE)
        self.assertEqual(info['seeds'], seeds)

    def test_size_gate(self):
        img = sitk.GetImageFromArray(self.test)

        blob_size = float(np.count_nonzero(self.test)) / \
                    reduce(lambda x, y: x*y, self.test.shape)

        (seeds, info) = sitkstrats.com_calc(img,
                                            blob_size - 0.001,
                                            0,
                                            self.ones_img)

        self.assertEqual(seeds, [])
        self.assertEqual(info['nseeds'], len(seeds))
        self.assertEqual(info['max_size'], blob_size - 0.001)

        (seeds, info) = sitkstrats.com_calc(img,
                                            self.MAX_SIZE,
                                            blob_size + 0.001,
                                            self.ones_img)

        self.assertEqual(seeds, [])
        self.assertEqual(info['nseeds'], len(seeds))
        self.assertEqual(info['min_size'], blob_size + 0.001)


    def test_mask(self):
        self.test[8:10, 6:10, 9:10] = 2
        img = sitk.GetImageFromArray(self.test)

        mask = np.zeros(self.test.shape)
        mask[5:10, 5:10, 5:10] = 1
        mask = sitk.GetImageFromArray(mask)

        (seeds, info) = sitkstrats.com_calc(img, self.MAX_SIZE,
                                            self.MIN_SIZE, mask)

        self.assertEqual(len(seeds), 1)
        self.assertEqual(len(seeds[0]), 3)
        self.assertEqual(type(seeds[0][0]), int)
        self.assertEqual(seeds[0], [9, 7, 8])

        (nmask_seeds, info) = sitkstrats.com_calc(img, 1e10, 0, self.ones_img)

        self.assertEqual(len(nmask_seeds), 2)
        self.assertEqual(len(nmask_seeds), info['nseeds'])
        self.assertEqual(nmask_seeds[1], seeds[0])
        self.assertEqual(nmask_seeds[0], [2, 2, 2])


if __name__ == '__main__':
    unittest.main()
