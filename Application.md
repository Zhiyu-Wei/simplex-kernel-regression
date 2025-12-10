# Application on WHOdataset

## Data description

During the COVID-19 pandemic, a widely discussed concern was whether educational attainment influenced people's hesitancy about vaccine side effects. To investigate this relationship, we utilized data from the \textit{2021 WHO Health Inequality Data Repository}, which provides the proportion of individuals at different education levels for each country. Since these values represent population shares, we treated them as compositional data, selecting three representative categories: 
at least secondary education, a four-year degree, and postgraduate education.
To further account for the role of institutional trust, we referenced the recent findings by \citet{Leblang2024}, who emphasized that trust in institutions significantly affects vaccination campaign outcomes. Guided by this insight, we incorporated two variables from the \textit{World Happiness Report (2021)}: index of institutional trust and perceptions of corruption.

After merging these sources, our final dataset consists of 72 countries, each characterized by compositional educational indicators and trust-related covariates. The response variable is defined as the proportion of unvaccinated individuals who reported concerns about side effects.

To explore potential nonlinear relationships between education compositions and side effect concerns, we visualize the data and model fits in Figure1. Based on the preceding simulation study, the LLD estimator demonstrates clear advantages when the data contain a substantial proportion of observations near the boundary of the simplex. As shown in Figure~\ref{fig:ObservedData}, the compositional component corresponding to Postgraduate education exhibits relatively low values (approximately 0.02–0.09), indicating that many observations lie close to the boundary. For this reason, we adopt the LLD method for the subsequent analysis.

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
```r
![Figure1](example/figures/concerned-plot.png)


