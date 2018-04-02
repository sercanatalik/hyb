import QuantLib as ql
from hyql.curves import rate
from hyql.utils import date
import pandas as pd

class Cashflow:
    """
     Creates a cashflow object, a payment a given date
    """
    def __init__(self, ccy, valDate, amount, sellBuyFlag, tradeDate):

        if str(type(valDate).__name__) == 'datetime':
            valDate = date.dt_to_QlDate(valDate)
        if str(type(tradeDate).__name__) == 'datetime':
            tradeDate = date.dt_to_QlDate(tradeDate)

        self.ccy = ccy
        self.amount = amount
        self.valDate = valDate
        self.tradeDate =tradeDate

        coupons = {'date': [],
                   'amount': [],
                   'ccy': []}

        coupons['date'].append(date.ql_to_DateTime(self.valDate))
        coupons['amount'].append(self.amount)
        coupons['ccy'].append(self.ccy)

        df = pd.DataFrame(coupons, columns=['date', 'amount','ccy'])
        df.index = df['date']
        del df['date']
        self.flows = df

    def npv(self, discountCurve):
        return self.amount * discountCurve.discount(self.valDate, True)

    



if __name__ == '__main__':
    today = ql.Date_todaysDate()
    ql.Settings.instance().evaluationDate = today
    cal = ql.UnitedStates()
    print(today)
    discountCurve = rate.build_flat_discounting_curve(0.1)
    c = Cashflow('USD', cal.advance(today, 1, ql.Years), 10e6, 'buy', today)
    print(c.npv(discountCurve))
    print(c.flows)