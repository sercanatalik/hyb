# coding=utf-8

from hyql.utils import date

import QuantLib as ql
import pandas as pd
import numpy as np


if __name__ == '__main__':

    today = ql.Date_todaysDate()

    date_grid = [today + ql.Period(i,ql.Months) for i in range(0,12*10)]

    date_grid = np.unique(np.sort(date_grid))
    time_grid = np.vectorize(lambda x: ql.ActualActual().yearFraction(today, x))(date_grid)
    dt = time_grid[1:] - time_grid[:-1]

    print(len(time_grid)*1500*2*29e-6)
    print(time_grid)

    seed = 1
    urng = ql.MersenneTwisterUniformRng(seed)
    usrg = ql.MersenneTwisterUniformRsg(len(time_grid) - 1, urng)
    generator = ql.InvCumulativeMersenneTwisterGaussianRsg(usrg)

    # Generate N paths
    N = 1500
    x = np.zeros((N, len(time_grid)))






