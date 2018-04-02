
import QuantLib as ql

from datetime import datetime

def ql_to_DateTime(d):
    return datetime(d.year(),d.month(),d.dayOfMonth())

def dt_to_QlDate(d):
    return ql.Date(d.day,d.month,d.year)

