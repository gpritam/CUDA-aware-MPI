## CUDA-aware-MPI

![CUDA-aware-MPI architecture](https://github.com/gpritam/CUDA-aware-MPI/blob/main/include/CPU-GPU-Architecture.png?raw=true)

## Dendritic growth

Dendritic growth is a phenomenon observed in fluid mechanics, specifically in the context of solidification processes. It refers to the formation of dendritic structures during the solidification of a liquid or molten material. When a liquid material undergoes solidification, it typically transitions from a disordered, liquid state to a more ordered, solid state. Dendritic growth occurs when the solidification process is characterized by the formation of tree-like or branch-like structures, known as dendrites. Dendrites are formed due to the anisotropic nature of the solidification process. As the liquid material cools and solidifies, it undergoes a phase transition, resulting in the formation of solid crystals. However, the growth of these crystals is not uniform in all directions. Instead, the crystals tend to grow preferentially along certain crystallographic directions, leading to the formation of dendritic structures. The growth of dendrites is influenced by various factors, including the temperature gradient, the concentration gradient, and the presence of impurities in the liquid material. These factors can affect the growth rate, morphology, and branching pattern of the dendrites.

Dendritic growth has significant implications in various fields, including materials science, metallurgy, and solidification processes. The morphology and microstructure of dendrites can impact the mechanical properties, thermal conductivity, and other characteristics of the solidified material. Therefore, understanding and controlling dendritic growth is crucial for optimizing the properties and performance of solidified materials. In fluid mechanics, the study of dendritic growth involves the analysis of the fluid flow and heat transfer during the solidification process. Computational fluid dynamics (CFD) simulations and mathematical models are often used to predict and analyze the dendritic growth patterns in different materials and processing conditions.

Here, we solve for two dependent variables namely the phase field parameter($\phi$) and the temperature($T$). The phase field parameter assumes the value 1 within the solid phase and 0 within the liquid phase. The governing equation for the phase field parameter is:

$$\tau \frac{\partial \phi}{\partial t} = \frac{\partial }{\partial y}\left(\epsilon \frac{\partial \epsilon}{\partial \theta} \frac{\partial \phi}{\partial x} \right) - \frac{\partial }{\partial x}\left(\epsilon \frac{\partial \epsilon}{\partial \theta} \frac{\partial \phi}{\partial y} \right) + \epsilon\nabla^2 \phi + \phi (1-\phi) \left(\phi + m - \frac{1}{2} \right)$$

The governing equation for the temperature becomes

$$\frac{\partial T}{\partial t} = \nabla^2 T + \kappa \frac{\partial \phi}{\partial t}$$

In the above equation, $\kappa$ is the dimensionless latent heat. The value of the outer normal vector at the solid-liquid interface is evaluated as 

$$\epsilon = \bar{\epsilon} \sigma,$$

where
$$\sigma = 1 + \delta \cos (\beta (\theta - \theta_0)),$$
and
$$\theta = \tan^{-1}\left(\frac{\partial \phi/\partial y}{\partial \phi/\partial x}\right).$$

The parameter $m$ is computed as below
$$m = \frac{\alpha}{\pi} \tan ^{-1} \left[ \gamma (T_{eq} - T)\right].$$

In our simulation, we assumed the following values for the constants: $\alpha = 0.9$, $\gamma=10.0$, $T_{eq}=1.0$, $\beta=6.0$, $\theta_0=0.2$, $\delta=0.02$, $\kappa=1.8$, and $\bar{\epsilon}=0.01$. We used the simplest explicit time integration scheme to solve these governing PDEs.


![Simulation of dendritic growth](https://github.com/gpritam/CUDA-aware-MPI/blob/main/Output/DendriticGrowth.gif?raw=true)

To build and run this code, invoke the following commands

```
make TARGET=DendriticGrowth_MPI.cpp
make run Nproc=1
```
In the above command, Nproc is the number of processors employed to solve the problem. The default value of this variable is 1. Upon successful completion of code execution, all the output data files in .tec format can be found in the Output directory. The output of the program can be visualized using Tecplot. At first, a shell script that is present in the Output directory is used to combine the results produced by different processors 
```
sh merge.sh
```
Next, images (in .png format) are generated in the same directory by running the script MovieDendriticGrowth.mcr in the Tecplot environment. Note that the |MFBD| variable in the .mcr script needs to be changed for the successful execution of the script. This particular variable denotes the path of the .tec data files. To generate the final movie file (in .gif format), invoke the following command

```
convert -resize 100% -delay 20 -loop 0 *.png DendriticGrowth.gif
```
The last command requires imagemagick to be installed. The result can be visualized in VisIt as well. To do so, we need to uncomment the following line in the DendriticGrowth.cpp file before the building process.
```
#define VISIT
```
Upon successful completion of code execution, all the output data files in .vtk format can be found in the Output directory. Open the Field.visit file in VisIt and play the animation.


## Spinoidal decomposition

Spinodal decomposition is a phenomenon that occurs in fluid mechanics, specifically in the study of phase separation in binary mixtures. It is a process where a homogeneous mixture spontaneously separates into two distinct phases due to the presence of thermodynamic instabilities. In a binary mixture, spinodal decomposition occurs when the system is in a metastable state, meaning it is not in equilibrium but is still stable for a certain period. This metastable state is characterized by a region in the phase diagram known as the spinodal region. The spinodal region is the region between the binodal curve and the spinodal curve. The binodal curve represents the boundary between the two-phase coexistence region and the single-phase region, while the spinodal curve represents the boundary between the metastable and unstable regions. When the system is in the spinodal region, small fluctuations in composition or density can grow exponentially, leading to the formation of distinct phases. These fluctuations are amplified by the negative second derivative of the free energy with respect to composition, which is a characteristic feature of the spinodal region.

In terms of applications, spinodal decomposition has significant implications in various fields, including materials science, polymer science, and colloid science. It plays a crucial role in the development of new materials with tailored properties, as well as in understanding the behavior of complex fluids.

We solve for the concentration field($C$) that is governed by the following equation
$$\frac{\partial C}{\partial t} = M\nabla^2 \left(\mu - \kappa \nabla^2 C\right),$$
where $\mu = 2\alpha C(1-C) (1-2C)$. In our simulation, we considered the following values for the constants: $\alpha=1$, $M=1$, and $\kappa = 0.5$.

![Simulation of spinoidal decomposition](https://github.com/gpritam/CUDA-aware-MPI/blob/main/Output/SpinoidalDecomposition.gif?raw=true)

To build and run this code, invoke the following commands

```
make TARGET=SpinoidalDecomposition_MPI.cpp
make run Nproc=1
```
In the above command, Nproc is the number of processors employed to solve the problem. The default value of this variable is 1. Upon successful completion of code execution, all the output data files in .tec format can be found in the Output directory. The output of the program can be visualized using Tecplot. At first, a shell script that is present in the Output directory is used to combine the results produced by different processors 
```
sh merge.sh
```
Next, images (in .png format) are generated in the same directory by running the script MovieSpinoidalDecomposition.mcr in the Tecplot environment. Note that the |MFBD| variable in the .mcr script needs to be changed for the successful execution of the script. This particular variable denotes the path of the .tec data files. To generate the final movie file (in .gif format), invoke the following command

```
convert -resize 100% -delay 10 -loop 0 *.png SpinoidalDecomposition.gif
```
The last command requires imagemagick to be installed. The result can be visualized in VisIt as well. To do so, we need to uncomment the following line in the SpinoidalDecomposition.cpp file before the building process.
```
#define VISIT
```
Upon successful completion of code execution, all the output data files in .vtk format can be found in the Output directory. Open the Field.visit file in VisIt and play the animation.

## References

1. S.B. Biner Programming Phase-Field Modeling. Springer
