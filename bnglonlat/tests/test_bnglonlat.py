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

        print('Testing inaccurate method')
        n = 100000
        xs = np.random.randint(0, 700000, size=n)
        ys = np.random.randint(0, 1250000, size=n)

        t0 = timeit.default_timer()
        aa, bb = bnglonlat.bnglonlat(xs, ys)
        t1 = timeit.default_timer() - t0

        t0 = timeit.default_timer()
        cc, dd = convertbng.util.convert_lonlat(xs, ys)
        t2 = timeit.default_timer() - t0

        print(f'Convertbng advantage: {round((t1 - t2) / t2 * 100)}%')

        e1 = np.abs(aa - cc)
        e2 = np.abs(bb - dd)
        ok = ~np.logical_or(np.isnan(e1), np.isnan(e2))
        e1 = e1[ok]
        e2 = e2[ok]
        print(f'Points returning nan: {n - len(e1)}')

        print(f'Max lon error: {np.max(e1)}')
        print(f'Max lat error: {np.max(e2)}')

        self.assertLess(np.max(e1), 5e-4)
        self.assertLess(np.max(e2), 2e-3)


if __name__ == '__main__':
    unittest.main()
