R code to generate optimal active portfolio allocations against a benchmark using a modified Grinold-Kahn (1999) framework.
It uses a risk-model based on historical factor-loadings/betas on excess-return factors.
In the current implemenetation, the alpha-signal is generated on a monthly basis and both the active and benchmark portfolios are rebalanced monthly.
Modify the signal_data.csv file in order to use your own alpha-signals.

**Current parameters**

*Benchmark Portfolio*

60% SPY
40% Barclays US Aggregate Bond Index

*Factor universe*

CAPM Market, Fama-French HML, SMB, UMD, CMA, RMW.

*Asset Universe*

Russell 1000 Growth and Value
Russell 2000 Growth and Value
Barclays US Aggregate Bond Index

*Transaction costs*

Provided in transaction_cost_data.csv

*Constraints*

Long-only, Max. Active risk < 15% annualized

*Beta significance-level*
p = 0.005 (two-tailed)
