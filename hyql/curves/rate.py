import QuantLib as ql


def build_flat_discounting_curve(flatRate,asofDate = ql.Settings.instance().evaluationDate):
    return ql.YieldTermStructureHandle(ql.FlatForward(asofDate, float(flatRate), ql.ActualActual()))





if __name__ == '__main__':
    today = ql.Date_todaysDate()
    ql.Settings.instance().evaluationDate = today
    print(today)

    discountCurve = rate.build_flat_discounting_curve(0.05)

    helpers = [ql.DepositRateHelper(ql.QuoteHandle(ql.SimpleQuote(rate / 100)),
                                    ql.Period(1, ql.Days), fixingDays,
                                    ql.TARGET(), ql.Following, False, ql.Actual360())
               for rate, fixingDays in [(0.04, 0), (0.04, 1), (0.04, 2)]]

    eonia = ql.Eonia()
    helpers += [ql.OISRateHelper(2, ql.Period(*tenor),
                                 ql.QuoteHandle(ql.SimpleQuote(rate / 100)), eonia)
                for rate, tenor in [(0.070, (1, ql.Weeks)), (0.069, (2, ql.Weeks)),
                                    (0.078, (3, ql.Weeks)), (0.074, (1, ql.Months))]]
    eonia_curve_c = ql.PiecewiseLogCubicDiscount(0, ql.TARGET(),
                                                 helpers, ql.Actual365Fixed())
    eonia_curve_c.enableExtrapolation()

    today = eonia_curve_c.referenceDate()
    end = today + ql.Period(2, ql.Years)
    dates = [ql.Date(serial) for serial in range(today.serialNumber(),
                                                 end.serialNumber() + 1)]
    rates_c = [eonia_curve_c.forwardRate(d, ql.TARGET().advance(d, 1, ql.Days),
                                         ql.Actual360(), ql.Simple).rate()
               for d in dates]

    print(rates_c)

