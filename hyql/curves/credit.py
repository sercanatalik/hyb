import QuantLib as ql


def generate_curveshape(spread, flat_curve_flag=False):
    m = [0, 100, 200, 500, 1000, 5000]
    ytwo = [0, 20, 60, 200, 500, 5000]
    yten = [0, 100, 250, 900, 1500, 5000]
    two = ql.LinearInterpolation(m, ytwo)
    ten = ql.LinearInterpolation(m, yten)
    if flat_curve_flag:
        return [spread, spread, spread]
    return [two(spread, True) / 10000, spread / 10000, ten(spread, True) / 10000]



def build_risky_Curve(discounting_curve, spreads,tenors=[ql.Period(10,ql.Years)],
                      asofDate=ql.Settings.evaluationDate(),cal = ql.UnitedStates,recoverRate = 0.40)
    instruments = [ql.SpreadCdsHelper(ql.QuoteHandle(ql.SimpleQuote(s)),
                                      tenor,
                                      0,
                                      cal,
                                      ql.Quarterly,
                                      ql.Following,
                                      ql.DateGeneration.TwentiethIMM,
                                      ql.Actual365Fixed(),
                                      recoverRate,
                                      discounting_curve)
                   for s,tenor in zip(spreads,tenors)]
    return ql.PiecewiseFlatHazardRate(asofDate,instruments,ql.Actual365Fixed())


