# Cyclostationary Spectral Proper Orthogonal Decomposition (CS-SPOD) in Matlab <a href="https://github.com/CyclostationarySPOD/CSSPOD/blob/main/LICENSE.md"> <img src="https://img.shields.io/badge/License-MIT-blue.svg" />


CSSPOD() is a Matlab implementation of the cyclostationary extension of spectral proper orthogonal decomposition (CS-SPOD) [1]. CS-SPOD is suitable for analyzing flows that exhibit periodic statistics. Examples of these cyclostationary flows include:		
	
  1. Flows with harmonic actuation (e.g. a periodically forced turbulent jet)	
  2. Turbines/rotors/VAWT/turbomachinery	
  3. Bluff body vortex shedding	
  4. Atmospheric (climate) data	
  5. Many others...
     
CS-SPOD determines the optimal decomposition of these cyclostationary flows in an analogous way that SPOD determines the optimal decomposition of statistically stationary flows. CS-SPOD becomes SPOD when analyzing statistically stationary flows. CS-SPOD is a multidimensional version of [3] using the method-of-snapshot technique for computational feasibility. 
</p>

### <p align="center"> <ins> To run the examples, download data from: https://caltech.box.com/v/CSSPODGithubData </ins> </p>

### Please cite the following publications if you use this code/method for your research:
<p align="center">
  <a href="https://arxiv.org/abs/2305.05628">
    Heidt, L. and Colonius, T., Journal of Fluid Mechanics (2024).
  </a>
</p>

```bibtex
@article{heidt2023spectral,   
     title={Spectral proper orthogonal decomposition of harmonically forced turbulent flows},   
     author={Heidt, Liam and Colonius, Tim},   
     journal={arXiv preprint arXiv:2305.05628},   
     year={2023}   
 }   
 ```

<p align="center">
  <a href="https://arc.aiaa.org/doi/10.2514/6.2023-3652">
    Heidt, L., Colonius, T., Nekkanti, A., Schmidt, O.T., Maia, I. and Jordan, P., AIAA AVIATION (2023).
  </a>
</p>
 
```bibtex
 @inproceedings{heidt2023cyclostationary,  
   title={Cyclostationary analysis of forced turbulent jets},    
   author={Heidt, Liam and Colonius, Tim and Nekkanti, Akhil and Schmidt, Oliver T and Maia, Igor and Jordan, Peter},  
   booktitle={AIAA AVIATION 2023 Forum},  
   pages={3652},  
   year={2023}  
 }
```

## Harmonically forced turbulent jet database
The large-eddy simulation data provided along with this example is a subset of the database of a harmonically forced (St_f=1.5) Mach 0.4 turbulent jet described in [2] and was calculated using the unstructured flow solver Charles developed at Cascade Technologies [4]. If you are using the database in your research or teaching, please include explicit mention of [2,4]. The test database consists of 25000 snapshots of the symmetric component (m=0) of a round turbulent jet. 

## References
[1] Heidt, L. and Colonius, T., 2023. Spectral proper orthogonal decomposition of harmonically forced turbulent flows. arXiv preprint arXiv:2305.05628. 

[2] Heidt, L., Colonius, T., Nekkanti, A., Schmidt, O.T., Maia, I. and Jordan, P., 2023. Cyclostationary analysis of forced turbulent jets. In AIAA AVIATION 2023 Forum (p. 3652).

[3] Kim, K.Y., North, G.R. and Huang, J., 1996. EOFs of one-dimensional cyclostationary time series: Computations, examples, and stochastic modeling. Journal of Atmospheric Sciences, 53(7), pp.1007-1017.

[4] Brès, Guillaume A., et al. "Importance of the nozzle-exit boundary-layer state in subsonic turbulent jets." Journal of Fluid Mechanics 851 (2018): 83-124.

## License
 
Copyright 2023.
This code is under the MIT license (see [LICENSE](LICENSE.md) file for full text).


