
# CVA calculation with QuantLib and Python

*Copyright (c) 2015 Matthias Groncki*


Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

    - Redistributions of source code must retain the above copyright notice,
    this list of conditions and the following disclaimer.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

This disclaimer is taken from the QuantLib license


```python
# import the used libraries
import numpy as np
import matplotlib.pyplot as plt
import QuantLib as ql
%matplotlib inline
```


```python
# Setting evaluation date
today = ql.Date(7,4,2015)
ql.Settings.instance().setEvaluationDate(today)
```


```python
# Setup Marketdata
rate = ql.SimpleQuote(0.03)
rate_handle = ql.QuoteHandle(rate)
dc = ql.Actual365Fixed()
yts = ql.FlatForward(today, rate_handle, dc)
yts.enableExtrapolation()
hyts = ql.RelinkableYieldTermStructureHandle(yts)
t0_curve = ql.YieldTermStructureHandle(yts)
euribor6m = ql.Euribor6M(hyts)
```


```python
# Setup a dummy portfolio with two Swaps
def makeSwap(start, maturity, nominal, fixedRate, index, typ=ql.VanillaSwap.Payer):
    """
    creates a plain vanilla swap with fixedLegTenor 1Y
    
    parameter:
        
        start (ql.Date) : Start Date
        
        maturity (ql.Period) : SwapTenor
        
        nominal (float) : Nominal
        
        fixedRate (float) : rate paid on fixed leg
        
        index (ql.IborIndex) : Index
        
    return: tuple(ql.Swap, list<Dates>) Swap and all fixing dates
    
        
    """
    end = ql.TARGET().advance(start, maturity)
    fixedLegTenor = ql.Period("1y")
    fixedLegBDC = ql.ModifiedFollowing
    fixedLegDC = ql.Thirty360(ql.Thirty360.BondBasis)
    spread = 0.0
    fixedSchedule = ql.Schedule(start,
                                end, 
                                fixedLegTenor, 
                                index.fixingCalendar(), 
                                fixedLegBDC,
                                fixedLegBDC, 
                                ql.DateGeneration.Backward,
                                False)
    floatSchedule = ql.Schedule(start,
                                end,
                                index.tenor(),
                                index.fixingCalendar(),
                                index.businessDayConvention(),
                                index.businessDayConvention(),
                                ql.DateGeneration.Backward,
                                False)
    swap = ql.VanillaSwap(typ, 
                          nominal,
                          fixedSchedule,
                          fixedRate,
                          fixedLegDC,
                          floatSchedule,
                          index,
                          spread,
                          index.dayCounter())
    return swap, [index.fixingDate(x) for x in floatSchedule][:-1]

portfolio = [makeSwap(today + ql.Period("2d"),
                      ql.Period("5Y"),
                      1e6,
                      0.03,
                      euribor6m),
             makeSwap(today + ql.Period("2d"),
                      ql.Period("4Y"),
                      5e5,
                      0.03,
                      euribor6m,
                      ql.VanillaSwap.Receiver),
            ]
```


```python
#%%timeit
# Setup pricing engine and calculate the npv
engine = ql.DiscountingSwapEngine(hyts)
for deal, fixingDates in portfolio:
    deal.setPricingEngine(engine)
    deal.NPV()
    #print(deal.NPV())
```


```python
# Stochastic Process 
```


```python
# Assume the model is already calibrated either historical or market implied
volas = [ql.QuoteHandle(ql.SimpleQuote(0.0075)),
         ql.QuoteHandle(ql.SimpleQuote(0.0075))]
meanRev = [ql.QuoteHandle(ql.SimpleQuote(0.02))]
model = ql.Gsr(t0_curve, [today+100], volas, meanRev, 16.)
```


```python
process = model.stateProcess()
```


```python
# Define evaluation grid
date_grid = [today + ql.Period(i,ql.Months) for i in range(0,12*6)]
for deal in portfolio:
    date_grid += deal[1]

date_grid = np.unique(np.sort(date_grid))
time_grid = np.vectorize(lambda x: ql.ActualActual().yearFraction(today, x))(date_grid)
dt = time_grid[1:] - time_grid[:-1]

print(len(time_grid)*1500*2*29e-6)
```

    6.786



```python
# Random number generator
seed = 1
urng = ql.MersenneTwisterUniformRng(seed)
usrg = ql.MersenneTwisterUniformRsg(len(time_grid)-1,urng)
generator = ql.InvCumulativeMersenneTwisterGaussianRsg(usrg)
```


```python
#%%timeit
# Generate N paths
N = 1500
x = np.zeros((N, len(time_grid)))
y = np.zeros((N, len(time_grid)))
pillars = np.array([0.0, 0.5, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
zero_bonds = np.zeros((N, len(time_grid), 12))

for j in range(12):
    zero_bonds[:, 0, j] = model.zerobond(pillars[j],
                                         0,
                                         0)
for n in range(0,N):
    dWs = generator.nextSequence().value()
    for i in range(1, len(time_grid)):
        t0 = time_grid[i-1]
        t1 = time_grid[i]
        x[n,i] = process.expectation(t0, 
                                     x[n,i-1], 
                                     dt[i-1]) + dWs[i-1] * process.stdDeviation(t0,
                                              x[n,i-1],
                                              dt[i-1])
        y[n,i] = (x[n,i] - process.expectation(0,0,t1)) / process.stdDeviation(0,0,t1)
        for j in range(12):
            zero_bonds[n, i, j] = model.zerobond(t1+pillars[j],
                                                 t1,
                                                 y[n, i])
                                                 
```


```python
# plot the paths
for i in range(0,N):
    plt.plot(time_grid, x[i,:])
```


![png](output_13_0.png)



```python
# generate the discount factors
discount_factors = np.vectorize(t0_curve.discount)(time_grid)
```


```python
#%%timeit
#Swap pricing under each scenario
npv_cube = np.zeros((N,len(date_grid), len(portfolio)))
for p in range(0,N):
    for t in range(0, len(date_grid)):
        date = date_grid[t]
        ql.Settings.instance().setEvaluationDate(date)
        ycDates = [date, 
                   date + ql.Period(6, ql.Months)] 
        ycDates += [date + ql.Period(i,ql.Years) for i in range(1,11)]
        yc = ql.DiscountCurve(ycDates, 
                              zero_bonds[p, t, :], 
                              ql.Actual365Fixed())
        yc.enableExtrapolation()
        hyts.linkTo(yc)
        if euribor6m.isValidFixingDate(date):
            fixing = euribor6m.fixing(date)
            euribor6m.addFixing(date, fixing)
        for i in range(len(portfolio)):
            npv_cube[p, t, i] = portfolio[i][0].NPV()
    ql.IndexManager.instance().clearHistories()
ql.Settings.instance().setEvaluationDate(today)
hyts.linkTo(yts)
```


```python
# Calculate the discounted npvs
discounted_cube = np.zeros(npv_cube.shape)
for i in range(npv_cube.shape[2]):
    discounted_cube[:,:,i] = npv_cube[:,:,i] * discount_factors
```


```python
# Calculate the portfolio npv by netting all NPV
portfolio_npv = np.sum(npv_cube,axis=2)
discounted_npv = np.sum(discounted_cube, axis=2)
```


```python
# Plot the first 30 NPV paths
n_0 = 0
n = 30
f, (ax1, ax2) = plt.subplots(2, 1, figsize=(12,10), sharey=True)
for i in range(n_0,n):
    ax1.plot(time_grid, portfolio_npv[i,:])
for i in range(n_0,n):
    ax2.plot(time_grid, discounted_npv[i,:])
ax1.set_xlabel("Time in years")
ax1.set_ylabel("NPV in time t Euros")
ax1.set_title("Simulated npv paths")
ax2.set_xlabel("Time in years")
ax2.set_ylabel("NPV in time 0 Euros")
ax2.set_title("Simulated discounted npv paths")
```




    <matplotlib.text.Text at 0x10ba4c810>




![png](output_18_1.png)



```python
# Calculate the exposure and discounted exposure
E = portfolio_npv.copy()
dE = discounted_npv.copy()
E[E<0] = 0
dE[dE<0] = 0
# Plot the first 30 exposure paths
n = 30
f, (ax1, ax2) = plt.subplots(2, 1, figsize=(12,10))
for i in range(0,n):
    ax1.plot(time_grid, E[i,:])
for i in range(0,n):
    ax2.plot(time_grid, dE[i,:])
ax1.set_xlabel("Time in years")
ax1.set_ylabel("Exposure")
ax1.set_ylim([-10000,70000])
ax1.set_title("Simulated exposure paths")
ax2.set_xlabel("Time in years")
ax2.set_ylabel("Discounted Exposure")
ax2.set_ylim([-10000,70000])
ax2.set_title("Simulated discounted exposure paths")
```




    <matplotlib.text.Text at 0x10bc10410>




![png](output_19_1.png)



```python
# Calculate the expected exposure
E = portfolio_npv.copy()
E[E<0]=0
EE = np.sum(E, axis=0)/N

```


```python
# Calculate the discounted expected exposure
dE = discounted_npv.copy()
dE[dE<0] = 0
dEE = np.sum(dE, axis=0)/N



```


```python
# plot the expected exposure path
n = 30
f, (ax1, ax2) = plt.subplots(2, 1, figsize=(8,10))
ax1.plot(time_grid, EE)
ax2.plot(time_grid, dEE)
ax1.set_xlabel("Time in years")
ax1.set_ylabel("Exposure")
ax1.set_title("Expected exposure")
ax2.set_xlabel("Time in years")
ax2.set_ylabel("Discounted Exposure")
ax2.set_title("Discounted expected exposure")

```




    <matplotlib.text.Text at 0x10c74ac90>




![png](output_22_1.png)



```python
# plot the expected exposure path
plt.figure(figsize=(7,5), dpi=300)
plt.plot(time_grid, dEE)
plt.xlabel("Time in years")
plt.ylabel("Discounting Expected Exposure")
plt.ylim([-2000,10000])
plt.title("Expected Exposure (netting set)")


```




    <matplotlib.text.Text at 0x103fca890>




![png](output_23_1.png)



```python
# Calculate the PFE curve (95% quantile)
PFE_curve = np.apply_along_axis(lambda x: np.sort(x)[0.95*N],0, E)

plt.figure(figsize=(7,5), dpi=300)
plt.plot(time_grid,PFE_curve)
plt.xlabel("Time in years")
plt.ylabel("PFE")
plt.ylim([-2000,35000])
plt.title("PFE (netting set)")
```




    <matplotlib.text.Text at 0x10405dfd0>




![png](output_24_1.png)



```python
# calculate the maximum pfe
MPFE = np.max(PFE_curve)
MPFE
```




    32034.275198581003




```python
# alternative pfe 95% quantile of the maxima of each exposure paths
PFE = np.sort(np.max(E,axis=1))[0.95*N]
PFE
```




    38683.65078838585




```python
# Setup Default Curve 
pd_dates =  [today + ql.Period(i, ql.Years) for i in range(11)]
hzrates = [0.02 * i for i in range(11)]
pd_curve = ql.HazardRateCurve(pd_dates,hzrates,ql.Actual365Fixed())
pd_curve.enableExtrapolation()
```


```python
# Plot curve
# Calculate default probs on grid *times*
times = np.linspace(0,30,100)
dp = np.vectorize(pd_curve.defaultProbability)(times)
sp = np.vectorize(pd_curve.survivalProbability)(times)
dd = np.vectorize(pd_curve.defaultDensity)(times)
hr = np.vectorize(pd_curve.hazardRate)(times)
f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(10,10))
ax1.plot(times, dp)
ax2.plot(times, sp)
ax3.plot(times, dd)
ax4.plot(times, hr)
ax1.set_xlabel("Time in years")
ax2.set_xlabel("Time in years")
ax3.set_xlabel("Time in years")
ax4.set_xlabel("Time in years")
ax1.set_ylabel("Probability")
ax2.set_ylabel("Probability")
ax3.set_ylabel("Density")
ax4.set_ylabel("HazardRate")
ax1.set_title("Default Probability")
ax2.set_title("Survival Probability")
ax3.set_title("Default density")
ax4.set_title("Harzard rate")


```




    <matplotlib.text.Text at 0x10b88b690>




![png](output_28_1.png)



```python
# Calculation of the default probs
defaultProb_vec = np.vectorize(pd_curve.defaultProbability)
dPD = defaultProb_vec(time_grid[:-1], time_grid[1:])
```


```python
# Calculation of the CVA
recovery = 0.4
CVA = (1-recovery) * np.sum(dEE[1:] * dPD)
CVA
```




    586.17763177944278


