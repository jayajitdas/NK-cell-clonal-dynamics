import pandas as pd
import numpy as np
import statsmodels.api as sm
import matplotlib.pyplot as plt
from scipy.stats import norm

# Example data
data = pd.read_csv('ly6c_cd27_data.csv')
x = data.Ly6c_percent.values
y = data.CD27_percent.values

# Add constant for intercept
X = sm.add_constant(x)
model = sm.OLS(y, X).fit()

# Predict values and get confidence intervals
x_pred = np.linspace(x.min(), x.max(), 100)
X_pred = sm.add_constant(x_pred)
preds = model.get_prediction(X_pred)
pred_summary = preds.summary_frame(alpha=0.05)  # 95% CI

from scipy.stats import pearsonr

r, p_value = pearsonr(x, y)
print(f"Correlation coefficient (r): {r}")
print(f"P-value: {p_value}")

def correlation_ci(r, n, alpha=0.05):
    # Fisher z-transformation
    z = np.arctanh(r)
    se = 1 / np.sqrt(n - 3)
    z_crit = norm.ppf(1 - alpha / 2)

    # CI in z-space
    z_low = z - z_crit * se
    z_high = z + z_crit * se

    # Transform back to r
    r_low = np.tanh(z_low)
    r_high = np.tanh(z_high)
    return r_low, r_high

n = len(x)
r_low, r_high = correlation_ci(r, n)

plt.scatter(x, y,c='k', label='Data')
plt.plot(x_pred, pred_summary['mean'], label='Regression line', color='red')
plt.fill_between(x_pred, pred_summary['mean_ci_lower'], pred_summary['mean_ci_upper'],
                 color='red', alpha=0.3, label='95% CI')
plt.legend()
plt.ylim(bottom=0)
plt.xlabel("%Ly6C+")
plt.ylabel("%CD27+")
plt.title(f"r = {r:.2f}, 95% CI: [{r_low:.2f}, {r_high:.2f}]")
plt.savefig('figure_ly6c_cd27.png')

