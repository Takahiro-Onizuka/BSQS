This repository provides R code implementing Gibbs sampler for Bayesian spatial quantile smoothing, as proposed by the following paper.

Onizuka, T., Hashimoto. S. and Sugasawa, S. (2022), Locally Adaptive Spatial Quantile Smoothing: Application to Monitoring Crime Density in Tokyo. *arXiv:2202.09534*.

The repository includes the following files.

* ```BSQS-function.R```: Implementation of Gibbs sampling for Bayesian spatial quantile smoothing methos: BQTF, SAR, and GP.

* ```BSQS-example-plot.R```: One-shot example of fitting Bayesian spatial quantile smoothing.

* ```TrueSignal-function.R```: The true signal functions such as two block structure and exponential function.

* ```NoiseDistribution-function.R```: The data-generating functions such as (I) Homogeneous, (II) Block heterogeneous, and (III) Smooth heterogeneous.
