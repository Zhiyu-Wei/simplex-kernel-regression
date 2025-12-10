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

<p align="center">
  <img src="figures/ori-plot.png" width="600">
  <br>
  <em>Figure 1: Observed data.</em>
</p>

The bandwidth parameter $h$ was selected using the leave-one-out cross-validation criterion,

$$
\hat{h} = \arg\min_{h \in (0,\infty)} \mathrm{CV}(h),
$$

which yielded an optimal value of approximately $\hat{h} = 0.0786$ when the Dirichlet kernel was used with tuning parameter $a = 0$. The choice of $a = 0$ is motivated by empirical performance: for this dataset, the estimator with $a = 0$ yields a smaller MSE compared to $a = 1$, and it also leads to more statistically decisive results in the subsequent model selection analysis.  

The figure below shows the graph of $CV(h)$ as a function of the bandwidth $h$.

<p align="center">
  <img src="figures/LLDbestCV.png" width="600">
  <br>
  <em>Figure 2: Plot of the leave-one-out cross-validation as a function of the bandwidth $h$.</em>
</p>

Then, using this selected bandwidth $\hat{h}$, we applied the LLD estimator on the simplex to visualize the fitted relationship and compared its performance with that of the linear model.

```r
#LLDK
composition_LL_kernel_smoother(data2,obj,0.005,h,original.data=data2,boundary.color="red")
#linear model
resolution<-0.005
Lmin <- min(data2[,1])
Rmin <- min(data2[,2])
Tmin <- min(data2[,3])

grid_points <- expand.grid(
  x = seq(0, 1, by = resolution),
  y = seq(0, 1, by = resolution)
)
grid_points$z <- 1 - grid_points$x - grid_points$y

grid_points <- subset(grid_points,x > Lmin & x < 1-Rmin-Tmin&y > Rmin & y < 1-Lmin-Tmin&z > Tmin & z < 1-Lmin-Rmin)
colnames(grid_points)<-c("x1","x2","x3")

model <- lm(y~x1+x2+x3,data=finaldata)

pred.lm <- predict(model,newdata=grid_points)
colnames(grid_points)<-colnames(data2)
draw_composition_plot(grid_points,pred.lm,boundary.color="red",original.data=data2,show_boundary = TRUE)
```

<p align="center">
  <img src="figures/WLS ori alp=0.png" width="350">
  <img src="figures/Linear.png" width="350">
  <br>
  <em>
    Figure 3: (Left) presents the fitted surface from Dirichlet kernel regression, capturing potential nonlinear
structures. (Right) shows the fitted surface from linear regression, which fails to capture the underlying structure
observed in the data.
  </em>
</p>

## Model selection
To formally evaluate whether nonlinear modeling improves the explanatory performance on real data, we performed pairwise model comparisons using the bootstrap likelihood ratio test (BLRT). This procedure tests whether including linear or kernel-smoothed compositional terms significantly enhances model fit relative to simpler alternatives.

**Table1** presents the BLRT results for three nested models. While adding linear compositional covariates **W1** vs **W2** does not yield a significant improvement, incorporating a nonparametric component **W3** leads to a significantly better fit when compared to both **W1** and **W2**.

<p align="center">
<strong>Table 1. BLRT results for model selection on real data.</strong> 
  <br>
<strong>The test evaluates whether more flexible models provide a significantly improved fit.</strong>
</p>

<table align="center">
  <tr>
    <th>Reduced Model</th>
    <th>Full Model</th>
    <th>p-value</th>
    <th>Reject H<sub>0</sub></th>
  </tr>
  <tr>
    <td><strong>W1</strong></td>
    <td><strong>W2</strong></td>
    <td>0.232</td>
    <td>No</td>
  </tr>
  <tr>
    <td><strong>W1</strong></td>
    <td><strong>W3</strong></td>
    <td>0.026</td>
    <td>Yes</td>
  </tr>
  <tr>
    <td><strong>W2</strong></td>
    <td><strong>W3</strong></td>
    <td>0.012</td>
    <td>Yes</td>
  </tr>
</table>

