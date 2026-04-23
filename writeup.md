---
title: "Density Time Series using Wasserstein Autoregression"
author: "Mason Griswold (Advised by Julia Fukuyama)"
date: "2026-04-16"
---

## Introduction

Traditionally, income distributions and longitudinal wealth dynamics are analyzed using aggregate statistics like the mean, median, or Gini coefficient.
While convenient, these scalar metrics obscure structural features of the underlying data, such as shifting multimodality, heavy right-skewness, and zero-inflation.
Capturing these dynamics requires moving beyond summary statistics and modeling the entire income distribution as a single functional object.

Functional Data Analysis (FDA) provides the tools to analyze sequences of continuous curves.
However, classical linear FDA methods fail when applied to probability density functions (PDFs) because they do not naturally enforce the requirements that densities remain non-negative and integrate to one.
Zhang et al. (2020) addressed this limitation by detailing specialized nonlinear transformations and geometric spaces that allow PDFs to be analyzed as valid data objects without violating these constraints.

Building on the time series frameworks discussed in Petersen et al. (2022), this paper benchmarks these methods against complex, real-world data.
While previous studies often evaluate density time series models on simulated or well-behaved datasets, we apply these frameworks to individual income data from the US Current Population Survey (CPS).
This data presents significant modeling challenges due to its extreme skewness and zero-inflation.

The primary objective of this study is to compare the predictive accuracy of competing Density Time Series methodologies: standard Functional Principal Component Analysis (FPCA), transformation-based approaches (LQD and Bayes spaces), and the geometric Wasserstein Autoregressive (WAR) model.
Although these frameworks support formal statistical inference, our scope is strictly restricted to evaluating their empirical forecasting performance.
By testing these models on weighted economic survey data, we determine which methodology is most effective when traditional assumptions of symmetry and normality fail.

## Functional Data Analysis

Conventional quantitative data typically takes the form of a scalar or a finite vector.
For instance, recording a subject's age, height, and weight yields an observation in $\mathbb{R}^3$.
However, in many applications, the underlying data-generating process is better represented as an entire continuous function.
For example, a subject's height recorded over several years represents a continuous growth curve.
While traditional time-series analysis could be applied to this data, it inherently assumes discrete time intervals.
Functional Data Analysis (FDA), by contrast, models the domain as a continuum, allowing for inference at arbitrary intermediate values without being constrained to discrete observation points (Ramsay and Silverman, 2005).

In practice, observed data is finite.
Because researchers rarely have the ability to measure a process continuously, a fundamental step in FDA is converting discrete data points into a functional curve defined over a continuum.
This is typically achieved using basis expansions, such as splines or Fourier series (Ramsay and Silverman, 2005).
A significant advantage of FDA is the ability to easily incorporate shape constraints.
For instance, human growth curves between ages 0 and 5 are strictly non-decreasing.
By enforcing this constraint during the smoothing process, FDA prevents the estimation of physically impossible functions that might otherwise arise from measurement error.

### Mathematical Framework and Construction

In standard FDA, we restrict our observed functions $f_i$ to a separable Hilbert space $\mathcal{H}$, most commonly the space of square-integrable functions, $L^2$ (Ramsay and Silverman, 2005).
Further restrictions depend on the specific analytical goals of the study.

Given a sample of $N$ underlying functions, each observation typically begins as a set of discrete data pairs $(x_{ij}, Y_{ij})$, where $j = 1, \dots, m_i$ represents the discrete measurements for the $i$-th subject.
From these discrete points, we construct the functional approximations $\hat{f}_i(t)$.
It is important to note that $\hat{f}_i(t)$ serves as an estimate of the true, unobservable continuous function $f_i(t)$, introducing some degree of error dependent on the measurement resolution and the chosen smoothing technique.

### Functional Principal Component Analysis

Once the functional approximations $\hat{f}_i(t)$ are constructed in a Hilbert space, exploratory analysis is often necessary to identify primary sources of variability and reduce the sample's dimensionality.
Functional Principal Component Analysis (FPCA) is the standard method for this task (Ramsay and Silverman, 2005).

FPCA relies on the Karhunen-Loève decomposition of a random functional element.
For a functional variable $Y(t)$ defined over a domain $\mathcal{T}$ with mean function $\nu(t) = \mathbb{E}[Y(t)]$, the decomposition is expressed as:

$$Y(t) = \nu(t) + \sum_{j=1}^{\infty} \xi_j \phi_j(t)$$

In this expansion:

- The $\phi_j(t)$ terms are orthogonal basis functions (eigenfunctions) that capture the most dominant directions of variability around the mean.
- The $\xi_j$ terms are the functional principal component scores, representing uncorrelated random variables with a mean of zero.

The scores are calculated by projecting the centered functions onto the basis functions:

$$\xi_j = \int_{\mathcal{T}} \{Y(t) - \nu(t)\}\phi_j(t)dt$$

By truncating this infinite sum to a finite number of $J$ components, we achieve dimension reduction.
Under the standard $L^2$ metric, this finite truncation is mathematically optimal in terms of the total variance explained (Ramsay and Silverman, 2005).

Because of its computational simplicity and highly interpretable outputs, FPCA is incredibly powerful for standard, unconstrained functional data.
However, when our functional data takes the form of objects with strict geometric rules, applying standard FPCA directly introduces significant mathematical contradictions.

### Functional Time Series and FPCA

In a conventional time series, observations are finite-dimensional scalars or vectors recorded at discrete time points.
A functional time series, however, sits at the intersection of time series analysis and functional data analysis: the overall sequence is still indexed by discrete time points ($t \in \mathbb{Z}$), but the observation recorded at each point is an entire continuous function residing in a function space.

This continuous nature presents a mathematical challenge for forecasting.
Because a continuous curve intrinsically possesses infinite dimensions, it cannot be reduced to a finite-dimensional space without utilizing basis expansions or specialized transformations.
Therefore, traditional forecasting models, such as Autoregressive (AR) or Vector Autoregressive (VAR) models, cannot be directly applied to the raw curves.

The common solution is to utilize FPCA.
By truncating this expansion to the first $J$ principal components, FPCA effectively compresses the infinite-dimensional curves into finite vectors of scores.
This unlocks the standard functional time series forecasting workflow:

1. Extract the basis functions and the corresponding FPC scores from the historical functional data.
2. Treat the sequence of extracted scores as a standard multivariate time series and fit a conventional model (such as a VAR model) to forecast the scores for future time steps.
3. Reconstruct the forecasted functional curves by multiplying the forecasted scores by the static basis functions.

While this workflow is highly effective for unconstrained data, it introduces significant issues when applied to probability density functions.

## Densities as Functional Data

Probability density functions (PDFs) represent a natural application for functional data analysis.
While parametric distributions can often be analyzed using standard multivariate techniques, non-parametric approaches are necessary when we wish to avoid imposing strict distributional assumptions on our data.
To construct our functional observations from raw sample points, we replace standard basis expansions (such as splines) with density estimation techniques, like Kernel Density Estimation (KDE), to generate our observed functions $\hat{f}_i$.

However, applying standard FDA techniques, such as FPCA, directly to these estimated densities introduces significant mathematical contradictions.
While well-behaved densities naturally reside in the Hilbert space $L^2$, probability densities are subject to strict non-linear constraints: they must be non-negative ($f \ge 0$) and integrate to one ($\int f = 1$).
Traditional linear FDA methods do not intrinsically respect these boundaries.
Because FPCA permits arbitrary linear combinations of principal component scores and basis functions, the resulting functional modes of variation are not guaranteed to remain bona fide probability density functions.
Consequently, reconstructing curves via standard FPCA can easily violate these fundamental probability constraints.

To overcome this limitation and properly model densities as data objects, modern functional data analysis relies on mapping these densities into representation spaces where linear methods become mathematically coherent (Petersen et al., 2022).
There are two primary schools of thought for achieving this: Transformation Approaches and Object-Oriented (Geometric) Approaches.

### Transformation Approaches

The transformation approach relies on applying an invertible mapping $\psi$ that projects a density function $f$ into a standard functional Hilbert space, such as $L^2$.
Once in this representation space, classical linear FDA techniques (like FPCA or linear regression) can be safely applied to the transformed curves.
The model outputs are then mapped back to the original density space $\mathcal{D}$ using the inverse transformation $\psi^{-1}$, which natively enforces the probability density constraints (Petersen and Muller, 2016).

Two prevalent transformations include:

- Log Quantile Density (LQD) Transformation: For a density $f$ with a corresponding quantile function $Q_f$, the LQD transformation is defined as $\psi_{LQD}(f)(t) = -\log\{f \circ Q_f(t)\}$ for $t \in [0,1]$ (Petersen and Muller, 2016).
  The representation space is $L^2[0,1]$, and the mathematical formulation of its inverse guarantees that the reconstructed function has a domain of $[0,1]$ and strictly integrates to 1.
  For densities without a common support (such as shifting income brackets over time), a modified LQD transformation can incorporate a location shift.
- Centered Log-Ratio (clr) Transformation: Rooted in compositional data analysis, this maps densities into a Bayes space, treating the probability density as carrying relative, rather than absolute, information (Kokoska et al., 2019).
  The clr transformation maps the densities to a subspace of $L^2$ containing functions that integrate to zero.

### Object-Oriented and Geometric Approaches (Wasserstein Space)

Alternatively, the object-oriented approach relies on the intrinsic geometric structure of the densities by utilizing non-Euclidean metrics that induce a manifold structure, such as the Fisher-Rao or Wasserstein metrics (Petersen et al., 2022).

When analyzing density time series, the 2-Wasserstein metric, an optimal transport distance that measures the cost of transforming one probability distribution into another, is highly effective.
Because the Wasserstein space $\mathcal{W}_2$ is nonlinear, linear time series models like autoregression cannot be applied directly to the density functions.
Instead, we can find a method to map the densities into a flat, linear space where traditional functional tools become mathematically valid.
This is achieved by operating within the tangent space.

The general framework for this geometric projection outlined in Zhang et al. 2020 proceeds as follows:

1. The Fréchet Mean $f_{\bigoplus}^{d_W}$: We first compute a central reference point for the sample of densities, known as the Wasserstein Fréchet mean.
   This minimizes the expected squared Wasserstein distance to all random densities in the sample and serves as our point of tangency.
   Note that Petersen et al. (2022) have found that the Fréchet mean is simple to calculate under the Wasserstein metric.
2. The Logarithmic Map: Observed densities $f_t$ are lifted into the tangent space at $f_\oplus$ using the logarithmic map, yielding $V_t = \text{Log}_{f_\oplus}(f_t)$.
3. Linear Modeling: Because the tangent space is a Hilbert space, it is perfectly flat and linear.
   Here, operations such as Tangent Space FPCA (Log-FPCA) or Wasserstein Autoregression (WAR) can be safely applied to the tangent vectors without violating probability constraints.
4. The Exponential Map: To transform a vector $V_*$ from this tangent space back into the Wasserstein space to become a density $f_*$, we use the exponential map.

## Forecasting Methods for Density Time Series

While the foundational mathematics of projecting densities into representation spaces is well-established, translating these geometries into predictive time series models has led to several distinct methodological approaches.
When observing a density time series $\{f_1, f_2, ..., f_n\}$, the goal is to forecast the next density, $f_{t+1}$, while preserving the probability constraints.

Recent literature has proposed a variety of models to tackle this, each with unique advantages depending on the data structure:

- Parametric Approaches (ST): Rather than modeling the full functional curve, this method assumes the data follows a specific parametric family, such as a skewed $t$-distribution.
  The parameters of the distribution are estimated at each time step, and traditional vector autoregression (VAR) is used to forecast the parameter vector.
  This is computationally efficient but highly restrictive if the true distributions exhibit complex, non-parametric behaviors.
- Dynamic Functional Principal Component Regression (HC): This method applies FPCA to the densities directly using a specific kernel and forecasts the resulting scores using VAR (Horta and Ziegelmann, 2018).
  Because this linear approach does not natively respect density constraints, any negative predicted values are artificially replaced by zero, and the reconstructed function is standardized to integrate to one.
- Log Quantile Density (LQD) Transformation: This relies on the LQD mapping to transport the densities into a Hilbert space where a multitude of functional data tools can be applied (Petersen and Muller, 2016).
  The inverse transformation is then applied to recover the forecasted density.
- Compositional Data Analysis (Bayes space): Rooted in the Bayes space geometry, this approach removes constraints by applying a centered log-ratio (clr) transformation (Kokoska et al. 2019).
  FPCA is applied to the transformed curves, and a time series model is fitted to the resulting coefficients.

While transformation-based methods like CoDa and LQD natively enforce the probability constraints of the predicted densities, they traditionally require the preliminary step of generating smoothed probability density functions (e.g., via Kernel Density Estimation) prior to transformation and analysis.

### The Wasserstein Autoregressive (WAR) Model

Rather than transforming the density functions into an unconstrained Hilbert space, the WAR model described by Zhang et al. (2020) operates directly on the tangent space of the Wasserstein manifold.

Once our densities are lifted into our linear tangent space, we may use an autoregressive model as such:

1. For a WAR($p$) model with scalar coefficients, the forecast is constructed via $\hat{V}_{t+1}(u) = \sum_{j=1}^p \hat{\phi}_j V_{t+1-j}(u)$, where $\hat{\phi}_j$ are the estimated autoregressive parameters.
2. The forecasted tangent vector $\hat{V}_{t+1}$ is projected back onto the density manifold using the exponential map to yield the final density forecast $\hat{f}_{t+1}$.

### Computational and Structural Advantages of WAR

The primary advantage of the WAR framework stems from the intrinsic properties of the 2-Wasserstein metric.
In one dimension, the Wasserstein distance between two measures is uniquely defined by the $L^2$ distance between their respective quantile functions (Zhang et al. 2020).
Because the entire geometric structure of $\mathcal{W}_2$ (including the Fréchet mean, tangent spaces, and logarithmic maps) can be expressed in terms of optimal transport maps between quantile functions, the WAR method effectively bypasses the need for explicit probability density functions during the modeling phase (Peterset et al. 2022).

This presents significant computational and methodological advantages over transformation-based methods like LQD or CoDa.
Specifically, the WAR model can be executed directly on cumulative distribution functions (CDFs) or empirical quantile grids.
By circumventing the need for Kernel Density Estimation (KDE), the WAR framework completely avoids the computational overhead of smoothing raw data points into functional curves.
Furthermore, this eliminates the risk of introducing estimation bias caused by subjective bandwidth parameter selection, allowing the autoregressive logic to interact directly with the empirical distribution of the data.
Structurally, the Wasserstein mean is also known to preserve the distinctive locations and heights of distributional modes far better than cross-sectional averaging or transformations.

### Application to Highly Skewed Survey Data

While previous research has benchmarked these density time series methods against well-behaved simulated datasets or financial returns, their comparative performance on highly non-normal, real-world data remains an open question.

This research aims to benchmark the predictive accuracy of the WAR method against transformation-based models (such as LQD and Bayes space) using longitudinal, individual-level income data derived from the Current Population Survey (CPS).
Income distributions introduce unique functional challenges: they are heavily right-skewed, naturally bounded at zero, and feature extreme outliers that heavily distort Euclidean metrics.
Furthermore, CPS data relies on complex survey weights to accurately reflect the true population structure, meaning the empirical quantiles and underlying densities are not uniform.

By analyzing how these competing geometries handle heavily skewed, weighted data, this study seeks to determine whether the computational and structural advantages of the quantile-driven Wasserstein framework yield superior forecasting accuracy when traditional assumptions of symmetry and normality break down.

## Methodology

To ensure a rigorous comparison between the KDE-based methods and the WAR framework, we require an objective and flexible bandwidth selection mechanism.
We utilize a $k$-nearest neighbors (KNN) adaptive bandwidth, which demonstrates robust performance for heavy-tailed distributions without compromising accuracy on standard distributions.

To determine the optimal parameter $k$, we employ a rolling-window forecasting approach.
For a candidate value of $k$, we estimate the probability density functions and construct a one-step-ahead FPCA forecast over a sequence of time points (e.g., iteratively forecasting $t=2000$ through $t=2024$).
Predictive accuracy is quantified by computing the 2-Wasserstein distance between the forecasted density and the true empirical distribution at each time step.
We evaluate multiple candidate values and select the $k$ that minimizes the average Wasserstein distance across the forecasting horizon.
To maintain consistency and isolate the effects of the respective transformations, this optimal $k$ is applied uniformly across all KDE-based models (Standard FPCA, LQD, and Bayes space).

In practice, functional data tools require continuous KDE functions to be discretized over a finite grid.
To achieve this, we evaluate the estimated densities over a pooled empirical quantile grid constructed from all observations across all time periods, utilizing probability increments of 0.001.
This approach guarantees that the resulting grid spans the full support of the data at all time points.
Furthermore, a quantile-based grid naturally allocates a higher density of evaluation points to regions with high probability mass while reducing resolution in the sparse tails.
This strategically optimizes the trade-off between functional fidelity and computational overhead, capturing critical structural details without the cost of a uniformly dense support grid.

To evaluate the accuracy of model predictions, we use the Wasserstein distance.
The Wasserstein distance between two 1-dimensional probability functions can be expressed as:

$$
W_2(f, g) = \left( \int_0^1 (Q_f(u) - Q_g(u))^2 \, du \right)^{1/2}
$$

$Q_f, Q_g$ are the quantile functions for $f$ and $g$ respectively (Zhang et al. 2020).

## Empirical Results: Synthetic Data

To benchmark the predictive accuracy of the models discussed, we first evaluated their performance on synthetic density time series.
For each scenario, we simulated $n = 1000$ independent observations at each discrete time point over a horizon of $T = 40$.
The data was generated from three distinct distributional families, Normal, Lognormal, and Uniform, to test the flexibility of each forecasting framework under varying geometric and boundary conditions.

Predictive accuracy was evaluated by computing the average 2-Wasserstein distance between the forecasted predictive densities and the true empirical quantiles of the holdout data.

### Normal Distribution

The normal distribution provides a baseline scenario with minimal structural complexity.
The mean value for each time point was generated via an autoregressive process:
$$\mu_t = 1.1\mu_{t-1} + \epsilon_t, \quad \epsilon_t \sim N(0, 0.10^2), \quad \mu_1 = 0$$
With the variance held constant, the distribution shifts symmetrically over time.

| Model                      | Average Wasserstein Distance |
| :------------------------- | :--------------------------- |
| FDA (Standard FPCA)        | 0.4870                       |
| Bayes (clr Transformation) | 0.5177                       |
| LQD Transformation         | 0.4939                       |
| WAR (Wasserstein AR)       | 0.5072                       |

![Quantile predictions of several Density Time Series methods on Normal data with an autoregressive mean.
Taken at the maximum time, 40](media/norm40.png)

Because the normal distribution is symmetric and lacks strict boundaries (defined on the entire real line), all methods perform comparably well.
Even unconstrained FPCA, which is less likely to violate probability constraints when the true density is situated far from zero does well here.
In scenarios where parametric assumptions are relaxed but the underlying data exhibits normal-like behavior, the WAR model represents a highly efficient choice, as it achieves equivalent accuracy without the overhead of applying Kernel Density Estimation (KDE).

### Lognormal Distribution

To introduce right-skewness and a strict lower bound at zero, we evolved a lognormal distribution using a similar autoregressive structure for the log-mean:
$$\mu_t = 0.98\mu_{t-1} + \epsilon_t, \quad \epsilon_t \sim N(0, 0.01^2), \quad \mu_1 = 10$$

| Model                      | Average Wasserstein Distance |
| :------------------------- | :--------------------------- |
| FDA (Standard FPCA)        | 42084.39                     |
| Bayes (clr Transformation) | 41882.08                     |
| LQD Transformation         | 42298.48                     |
| WAR (Wasserstein AR)       | 42077.70                     |

![Quantile predictions of several Density Time Series methods on Log-Normal data with an autoregressive mean.
Taken at the maximum time, 40](media/log_norm40.png)

Overall, the transformation-based models (LQD and Bayes) successfully capture the right-skewed geometry of the data.
While the standard FPCA approach tracks the distribution reasonably well in early periods, it becomes increasingly unstable over time.
This instability occurs because the support of the underlying distribution shifts significantly as the series progresses.
Since FPCA constructs forecasts strictly through linear combinations of historically observed basis functions, it intrinsically struggles to accurately estimate future probability densities that reside on novel or shifted domains.
The WAR model, by contrast, correctly captures the functional shape of the distribution, but contains a significant location shift.
This location shift could likely be solved by introducing differencing to the WAR model.

### Uniform Distribution

Finally, to test how the models handle strict, hard boundaries, we simulated a uniform distribution where the bounds $a_t$ and $b_t$ independently evolve over time:
$$a_t = 1.05a_{t-1} + \epsilon_{a,t}, \quad b_t = 1.05b_{t-1} + \epsilon_{b,t}$$
Where $a_1 = 0$, $b_1 = 1$, and $\epsilon_{a,t}, \epsilon_{b,t} \sim N(0, 0.10^2)$.

| Model                      | Average Wasserstein Distance |
| :------------------------- | :--------------------------- |
| FDA (Standard FPCA)        | 1.4586                       |
| Bayes (clr Transformation) | 1.8044                       |
| LQD Transformation         | 1.1432                       |
| WAR (Wasserstein AR)       | 1.1992                       |

![Quantile predictions of several Density Time Series methods on Uniform data with an autoregressive lower and upper bounds.
Taken at the maximum time, 40](media/unif40.png)

The results here highlight the methodological differences between the frameworks.
Standard FPCA and the Bayes space approach perform remarkably poorly.
The Bayes space transformation introduces extreme visual distortions, creating artificial "humps" in the early quantiles.

Part of this failure is directly attributable to the preliminary density estimation step.
Transformation models require raw data to be smoothed into PDFs via KDE before applying the transformation.
When a Gaussian kernel is applied to uniform data, it forces smooth, asymptotic tails beyond the strict $[a, b]$ boundaries.
It would be optimal to use a different choice of kernel for densities with strict cutoffs.
A solution to this issue would likely make the LQD transformation a clear winner.
It exhibits higher predictive accuracy, but currently suffers from smoothly decaying tails.

We are unsure of why the Bayes space method does poorly here.
The method does quite well for some time points, but quickly deviates to include the aforementioned "hump".
This instability seems to be isolated to the Uniform distribution.

The WAR method, by operating directly on the empirical cumulative distribution functions (CDFs) rather than smoothed PDFs, respects the hard boundaries of the uniform distribution without requiring the subjective tuning of a kernel bandwidth.
The WAR also exhibits good predictive accuracy with a Wasserstein distance within 5% of the LQD method.

## Real Data: US Individual Incomes

Evaluating predictive models exclusively on synthetic datasets may not accurately reflect their performance in real-world applications.
To assess their practical utility, we apply the forecasting frameworks to US individual income data spanning from 1998 to 2024.
This data is sourced from the US Census Bureau's Current Population Survey (CPS), a comprehensive annual survey of American households.

Our analysis focuses on two specific variables from this survey: `ptotval` and `a_ernlwt`.
The `ptotval` variable records total individual income in nominal US dollars.
To ensure comparability across the temporal domain, these values were manually scaled using Consumer Price Index (CPI) data, converting them into inflation-adjusted real dollars.
The `a_ernlwt` variable provides the corresponding survey weights for the earnings questions.
Because the CPS utilizes a complex sampling design, incorporating these survey weights to construct weighted empirical distributions is essential to accurately reflect the true population structure.

Beyond the integration of these survey weights, the forecasting methodology aligns identically with the structure utilized in the synthetic data analysis.

| Model                      | Average Wasserstein Distance |
| :------------------------- | :--------------------------- |
| FDA (Standard FPCA)        | 149754.1                     |
| Bayes (clr Transformation) | 137338.2                     |
| LQD Transformation         | 135701.0                     |
| WAR (Wasserstein AR)       | 150603.0                     |

Quantitatively, the transformation-based models (LQD and Bayes) achieved slightly lower average Wasserstein distances than the WAR model.
However, visual inspection reveals that they provide a noticeably poorer structural match to the empirical data.
In particular, the LQD transformation exhibits an extreme deviation from the observed quantiles at probability values below 0.5.
The method struggles to capture the zero-inflated nature of the income distribution.

Conversely, by bypassing the KDE step and operating directly on the weighted empirical quantiles, the WAR method visually outperforms its counterparts, successfully preserving the zero-inflated geometry of the income data.

![Quantile predictions of several Density Time Series methods on CPS Income data.
Taken at the maximum time, 2024](media/incomes2024.png)

## Conclusion

This study compared the predictive accuracy of density time series forecasting models, contrasting transformation-based approaches (LQD and Bayes space) with the Wasserstein Autoregressive (WAR) framework.
While conventional Functional Principal Component Analysis (FPCA) inherently violates the non-negative and integrability constraints of probability density functions, modern approaches navigate these boundaries by projecting data into specialized representation spaces.

Our empirical results highlight a critical bottleneck in transformation-based methodologies: their reliance on a preliminary Kernel Density Estimation (KDE) step.
As demonstrated by the synthetic uniform distribution and the CPS individual income data, KDE adds complexity but generally performs on par with WAR methods.
By operating directly on weighted empirical quantiles and bypassing continuous density approximation, the WAR model avoids KDE-induced boundary effects.
Even when LQD and Bayes approaches yielded marginally lower average error metrics, the WAR framework better preserved key structural features, such as heavy tails and zero-inflated peaks, in highly skewed socioeconomic data.

Despite these advantages, the standard WAR framework has limitations.
Our analysis of the synthetic lognormal distribution showed that the model can exhibit location shifts when applied to non-stationary data with dynamically evolving means.
Future research should explore differenced Wasserstein autoregressive models to account for these non-stationary level shifts.
Additionally, expanding this comparative framework to include Fully Functional WAR (FF-WAR) models could refine our understanding of optimal transport geometries in time series forecasting.

For researchers analyzing strictly bounded, heavily skewed, or complex weighted data, the quantile-driven WAR model offers a structurally faithful and computationally efficient forecasting solution.

## References

Delicado, P. ‘Dimensionality Reduction When Data Are Density Functions’. Computational Statistics & Data Analysis, vol. 55, no. 1, 2011, pp. 401–420, <https://doi.org/10.1016/j.csda.2010.05.008>.

Horta, Eduardo, and Flavio Ziegelmann. ‘Dynamics of Financial Returns Densities: A Functional Approach Applied to the Bovespa Intraday Index’. International Journal of Forecasting, vol. 34, no. 1, 2018, pp. 75–88, <https://doi.org/10.1016/j.ijforecast.2017.08.001>.

Kokoszka, Piotr, et al. ‘Forecasting of Density Functions with an Application to Cross-Sectional and Intraday Returns’. International Journal of Forecasting, vol. 35, no. 4, 2019, pp. 1304–1317, <https://doi.org/10.1016/j.ijforecast.2019.05.007>.

Petersen, Alexander, and Hans-Georg Müller. ‘Functional Data Analysis for Density Functions by Transformation to a Hilbert Space’. The Annals of Statistics, vol. 44, no. 1, Institute of Mathematical Statistics, 2016, pp. 183–218, <https://doi.org/10.1214/15-AOS1363>.

Petersen, Alexander, et al. ‘Modeling Probability Density Functions as Data Objects’. Econometrics and Statistics, vol. 21, 2022, pp. 159–178, <https://doi.org/10.1016/j.ecosta.2021.04.004>.

Ramsay, James O., and Bernard W. Silverman. Functional Data Analysis. Springer, 1997.

Zhang, Chao, et al. ‘Wasserstein Autoregressive Models for Density Time Series’. arXiv [Stat.ME], 2020, arxiv.org/abs/2006.12640. arXiv.
