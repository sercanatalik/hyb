import QuantLib as ql
from hyql.curves import rate
from hyql.utils import date
import pandas as pd
from hyql.products import cds


class MultiAssetBook:
    def __init__(self):
        self.cdsPorfolio = []
        self.fxPortfolio = []

    def add_cds(self,cdsTrade):
        self.cdsPorfolio.append(cdsTrade)

    def add_fx(self,fxTrade):
        self.fxPortfolio.append(fxTrade)

    def set_cdsEngine(self,engine):
        for c in self.cdsPorfolio:
            c.set_pricingengine(engine)






if __name__ == '__main__':
    print('testing portfolio')

    issuer = 'TURKEY'
    startDate = ql.Date_todaysDate()
    endDate = ql.Date_todaysDate()+360

    cds1 = cds.CDS(issuer, startDate, endDate, 100/10000, 10e6,'Buy')

    bok = MultiAssetBook()
    bok.add_cds(cds1)

    print(bok.cdsPorfolio)

    print(cds1)