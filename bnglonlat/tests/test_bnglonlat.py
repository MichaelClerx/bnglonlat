#!/usr/bin/env python3
#
# Tests for bnglonlat
#
import unittest
import timeit

import convertbng.util
import numpy as np

import bnglonlat


class TestBngLonlat(unittest.TestCase):
    """ Tests for ``bnglonlat``. """

    def test_random(self):
        # Compares random points against convertbng

        n = 100000
        xs = np.random.randint(0, 700000, size=n)
        ys = np.random.randint(0, 1300000, size=n)

        t0 = timeit.default_timer()
        aa, bb = bnglonlat.bnglonlat(xs, ys)
        t1 = timeit.default_timer() - t0

        t0 = timeit.default_timer()
        cc, dd = convertbng.util.convert_lonlat(xs, ys)
        t2 = timeit.default_timer() - t0

        print(f'Convertbng advantage: {round((t1 - t2) / t2 * 100)}%')

        e1 = np.abs(aa - cc)
        e2 = np.abs(bb - dd)
        e1 = e1[~np.isnan(e1)]
        e2 = e2[~np.isnan(e2)]

        print(f'Max lon error: {np.max(e1)}')
        print(f'Max lat error: {np.max(e2)}')

        self.assertLess(np.max(e1), 5e-4)
        self.assertLess(np.max(e2), 2e-3)


if __name__ == '__main__':
    unittest.main()
