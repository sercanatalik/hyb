# coding=utf-8

from hyql.utils import date

import QuantLib as ql
import pandas as pd
import numpy as np

if __name__ == '__main__':
    today = ql.Date_todaysDate()

    # Setup Marketdata
    rate = ql.SimpleQuote(0.03)
    rate_handle = ql.QuoteHandle(rate)
    dc = ql.Actual365Fixed()
    yts = ql.FlatForward(today, rate_handle, dc)
    yts.enableExtrapolation()
    hyts = ql.RelinkableYieldTermStructureHandle(yts)
    t0_curve = ql.YieldTermStructureHandle(yts)
    euribor6m = ql.Euribor6M(hyts)
    


    date_grid = [today + ql.Period(i, ql.Months) for i in range(0, 12 * 10)]

    date_grid = np.unique(np.sort(date_grid))
    time_grid = np.vectorize(lambda x: ql.ActualActual().yearFraction(today, x))(date_grid)
    dt = time_grid[1:] - time_grid[:-1]

    print(len(time_grid) * 1500 * 2 * 29e-6)
    print(time_grid)

    seed = 1
    urng = ql.MersenneTwisterUniformRng(seed)
    usrg = ql.MersenneTwisterUniformRsg(len(time_grid) - 1, urng)
    generator = ql.InvCumulativeMersenneTwisterGaussianRsg(usrg)

    # Generate N paths
    N = 10
    x = np.zeros((N, len(time_grid)))

    print(generator.dimension())

    npv_cube = np.zeros((N, len(date_grid), 1)) #portfolio=1

    for p in range(0, N):
        for t in range(0, len(date_grid)):
            date = date_grid[t]
            print(date)
            ql.Settings.instance().setEvaluationDate(date)
            ycDates = [date,
                       date + ql.Period(6, ql.Months)]
            ycDates += [date + ql.Period(i, ql.Years) for i in range(1, 11)]

            print(ycDates)
         #   yc = ql.DiscountCurve(ycDates,
          #                        zero_bonds[p, t, :],
           #                       ql.Actual365Fixed())
          #  yc.enableExtrapolation()
          #  hyts.linkTo(yc)

           # for i in range(len(portfolio)):
           #     npv_cube[p, t, i] = portfolio[i][0].NPV()
        #ql.IndexManager.instance().clearHistories()
    ql.Settings.instance().setEvaluationDate(today)
    #hyts.linkTo(yts)