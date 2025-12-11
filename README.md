## Figures

This section contains the figures generated for the paper **Kernel Regression with Simplex Predictors**.  
All visualizations were produced using the custom ternary plotting function implemented in this repository.  
The function `draw_composition_plot()` maps compositional predictors on the 2-simplex to barycentric coordinates and overlays kernel-smoothed responses, allowing clear visualization of nonlinear regression surfaces and boundary behavior.

The figures below illustrate the regression functions used in the simulation study, evaluated across a dense simplex grid.  
Each plot corresponds to one of the true functions \( m_1, m_2, m_3, m_4 \), and demonstrates how the response varies over the simplex domain.

### **Visualization Function**
- **`draw_composition_plot()`**  
  - Converts 3-component compositions to 2-dimensional barycentric coordinates  
  - Draws simplex boundaries and optional grid structure  
  - Supports continuous color gradients for response visualization  
  - Designed for kernel regression diagnostics and shape exploration  
  - Outputs publication-quality ternary plots via `CairoPNG`

### **Figures (generated from the code `plot_m1-m4.R`):**

#### **m1: Linear combination of simplex components**

<p align="center">
  <img src="m1.png" width="600">
  <br>
</p>


#### **m2: Logarithmic surface**

<p align="center">
  <img src="m2.png" width="600">
  <br>
</p>

#### **m3: Interaction surface \(x_1 x_2\)**
<p align="center">
  <img src="m3.png" width="600">
  <br>
</p>

#### **m4: Nonlinear oscillating structure**

<p align="center">
  <img src="m4.png" width="600">
  <br>
</p>
