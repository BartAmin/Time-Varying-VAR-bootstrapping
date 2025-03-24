# Time-Varying VAR Bootstrap NHST

Understanding whether the underlying stochastic process of a time series is (weakly) stationary is crucial before applying traditional time series models. If the process is stationary, models such as ARIMA, VAR, or GARCH can be applied reliably. However, applying these models to data from non-stationary processes may lead to biased or misleading parameter estimates. In such cases, time-varying models—such as Time-Varying Vector Autoregression (tv-VAR)—offer a more appropriate alternative.
One way to evaluate whether a tv-VAR model is warranted is through a parametric bootstrap hypothesis test. This involves simulating multiple stationary time series and fitting tv-VAR models to these synthetic datasets to estimate the expected variance of parameters under stationarity. Next, a tv-VAR model is fitted to the original dataset, and the variance of its parameters is compared to the bootstrapped null distribution. This comparison yields p-values that help determine whether the observed parameter variability is greater than expected under stationarity.
If the p-value is below a chosen significance level (e.g., α = 0.05), it suggests that the underlying process is likely non-stationary and that a tv-VAR model may be more appropriate. This script evaluates the statistical properties—such as Type I error rate and power—of this approach to assess its reliability.

## The script:
- Simulates both stationary and non-stationary multivariate time series data
- Fits a VAR(1) model to extract time-invariant parameters
- Simulates new datasets using these estimated parameters
- Applies the `tvvarGAM` method to the synthetic datasets
- Computes the variability (standard deviation) of estimated parameters across bootstrap iterations
- Calculates p-values by comparing the variability in the original dataset to the bootstrapped variability
- Evaluates Type I error control using ECDF plots
- Assesses classification performance using ROC curves and AUC

### Dependencies

This project requires the following R packages:

- [`mlVAR`](https://cran.r-project.org/package=mlVAR)  
- [`mgm`](https://cran.r-project.org/package=mgm)  
- [`devtools`](https://cran.r-project.org/package=devtools) 
- [`tvvarGAM`](https://github.com/LauraBringmann/tvvarGAM) 
- [`pROC`](https://cran.r-project.org/package=pROC) 
- [`ggplot2`](https://cran.r-project.org/package=ggplot2)

To install `tvvarGAM` from GitHub:

```r
install.packages("remotes")
remotes::install_github("LauraBringmann/tvvarGAM")
```

### Acknowledgement
- For a full description of how the Time-Varying VAR with Generalized Additive Models works please consult: https://github.com/LauraBringmann/tvvarGAM
- For a full description of how the mixed-VAR works please consult: https://github.com/jmbh/mgm
