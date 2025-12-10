# Application on WHOdataset

## Data description

During the COVID-19 pandemic, a widely discussed question was whether educational attainment influenced people's hesitancy about vaccine side effects. To investigate this relationship, we used data from the **2021 WHO Health Inequality Data Repository**, which reports the proportion of individuals at different education levels for each country. Since these values represent population shares, we treated them as compositional data and focused on three representative categories: at least secondary education, a four-year degree, and postgraduate education.

To incorporate the potential role of institutional trust, we also included two variables from the **World Happiness Report (2021)**: institutional trust and perceptions of corruption. Earlier studies have suggested that institutional trust plays an important role in shaping public responses to vaccination campaigns, and these variables help capture that dimension.

After merging these sources, the final dataset contains 72 countries, each described by compositional educational indicators and trust-related covariates. The response variable is the proportion of unvaccinated individuals who reported concerns about side effects.

To explore possible nonlinear relationships between education compositions and concerns about side effects, we visualize the data and model fits in the figure below. Prior simulation studies showed that the local linear Dirichlet (LLD) estimator performs particularly well when compositional observations lie near the boundary of the simplex. As seen in the observed data, the proportion corresponding to postgraduate education is relatively small (approximately 0.02–0.09), indicating that many observations lie close to the boundary. Motivated by this feature, we adopt the LLD method for the subsequent analysis.

```r
library(doParallel)
library(ggplot2)
library(viridis)
source(""https://raw.githubusercontent.com/Zhiyu-Wei/simplex-kernel-regression/R/application.R")
url <- "https://raw.githubusercontent.com/Zhiyu-Wei/simplex-kernel-regression/example/WHOdataset.csv"
ds <- read.csv(url)
obj=ds[,3]
data1=ds[,c(10,19)]
data2=as.matrix(ds[,c(6,7,8)])
finaldata <- cbind(obj,data2,data1)
colnames(finaldata)=c("y","x1","x2","x3","z1","z2")
n <- nrow(finaldata)
draw_composition_plot(data2,obj,boundary.color="red",original.data=data2)
```
<div align="center">
<img src="figures/ori-plot.png" width="600">
*Figure 1: Observed data.*
</div>

<img src="figures/observed_education_simplex.png" width="600">


The bandwidth parameter $h$ was selected using the leave-one-out cross-validation criterion,

$$
\hat{h} = \arg\min_{h \in (0,\infty)} \mathrm{CV}(h),
$$

which yielded an optimal value of approximately $\hat{h} = 0.0786$ when the Dirichlet kernel was used with tuning parameter $a = 0$. The choice of $a = 0$ is motivated by empirical performance: for this dataset, the estimator with $a = 0$ yields a smaller MSE compared to $a = 1$, and it also leads to more statistically decisive results in the subsequent model selection analysis.  

The figure below shows the graph of $CV(h)$ as a function of the bandwidth $h$.

![Figure2](figures/LLDbestCV.png)

Then, using this selected bandwidth $\hat{h}$, we applied the LLD estimator on the simplex to visualize the fitted relationship and compared its performance with that of the linear model.


