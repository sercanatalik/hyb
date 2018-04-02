import QuantLib as ql
import pandas as pd
from hyql.utils import date


def cds_singlename(startDate, endDate, spread, Notional, payReceiveFlag):
    if str(type(startDate).__name__) == 'datetime':
        startDate = date.dt_to_QlDate(startDate)
    if str(type(endDate).__name__) == 'datetime':
        endDate = date.dt_to_QlDate(endDate)
    cal = ql.UnitedStates()
    if payReceiveFlag == 'Buy':
        flag = ql.Protection.Buyer
    if payReceiveFlag == 'Sell':
        flag = ql.Protection.Seller

    schedule = ql.Schedule(startDate, endDate, ql.Period(ql.Quarterly),
                           cal, ql.Following, ql.Unadjusted, ql.DateGeneration.TwentiethIMM, False)
    return ql.CreditDefaultSwap(flag, Notional, spread, schedule, ql.Following, ql.Actual365Fixed())


def cds_coupons(cds):
    coupons = {'date': [],
               'amount': []}
    c = cds.coupons()
    for i in range(len(c)):
        coupons['date'].append(date.ql_to_DateTime(c[i].date()))
        coupons['amount'].append(c[i].amount())
    df = pd.DataFrame(coupons, columns=['date', 'amount'])
    df.index = df['date']
    del df['date']
    return df


class CDS:
    def __init__(self, issuer, startDate, endDate, spread, notional, payReceiveFlag):
        if str(type(startDate).__name__) == 'datetime':
            startDate = date.dt_to_QlDate(startDate)
        if str(type(    endDate).__name__) == 'datetime':
            endDate = date.dt_to_QlDate(endDate)
        self.issuer = issuer
        self.startDate = startDate
        self.endDate = endDate
        self.notional = notional
        self.payReceiveFlag = payReceiveFlag
        self.ql_cds = cds_singlename(startDate, endDate, spread, notional, payReceiveFlag)
        self.premium_paid = 0
        self.pricingEngine = None
        self.side = self.ql_cds.side()
        self.lastCalculatedNpv = 0
        self.lastCalculatedPnl = 0

    def maturity(self):
        c = self.ql_cds.coupons()
        return c[-1].date()

    def accrualStart(self):
        c = self.ql_cds.coupons()
        return c[0].date()

    def cs01(self):
        return -1 * self.ql_cds.couponLegBPS()

    def forward(self):
        return self.ql_cds.fairSpread()

    def npv(self):
        if self.pricingEngine is None:
            return 'error : set pricing engine for cds  ' + self.issuer
        npv = self.ql_cds.NPV()
        self.lastCalculatedNpv = npv
        return npv

    def set_pricingengine(self, engine):
        self.pricingEngine = engine
        self.ql_cds.setPricingEngine(self.pricingEngine)

    def get_pricingengine(self):
        return self.pricingEngine

    def forward_spread(self):
        if self.pricingEngine is None:
            return 'error : set pricing engine for cds  ' + self.issuer
        return self.ql_cds.fairSpread()

    def get_coupons(self):
        return cds_coupons(self.ql_cds)

    def get_past_coupons(self, asofDate):
        c = self.get_coupons()
        return c[:date.ql_to_DateTime(asofDate)]

    def get_future_coupons(self, asofDate):
        c = self.get_coupons()
        return c[date.ql_to_DateTime(asofDate):]

    def set_premium_paid(self, premium):
        self.premium_paid = premium

    def get_past_carry(self, asofDate=ql.Settings.evaluationDate):
        amt = self.get_past_coupons(asofDate).sum()['amount']
        if self.ql_cds.side() == 0:
            amt = -1 * amt
        return amt

    def get_pnl(self,asofDate=ql.Settings.evaluationDate):
        pnl = self.npv()+self.premium_paid+self.get_past_carry(asofDate)
        self.lastCalculatedPnl = pnl
        return pnl

