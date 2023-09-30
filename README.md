# HyRes_GPU model tutorial

## Shortly
HyRes model is a **Hy**brid **Res**olution coarse-grained protein model developed since 2017, where the backbone is atomistic-style description but coarse-grained side chains. For a full and complete description of the HyRes model, please move to the related papers as follows:  
1. HyRes: [https://doi.org/10.1039/C7CP06736D](https://doi.org/10.1039/C7CP06736D)
2. HyRes-II: [https://doi.org/10.1021/acs.jcim.2c00974](https://doi.org/10.1021/acs.jcim.2c00974)
3. HyRes-GPU: [https://doi.org/10.1101/2023.08.22.554378](https://doi.org/10.1101/2023.08.22.554378)  

For running a HyRes simulation, two steps are necessary:
1. convert atomistic model to HyRes model ([at2hyres](https://github.com/wayuer19/HyRes_GPU/tree/main/at2hyres))
2. describe HyRes force field in OpenMM and run simulation ([run_OpenMM](https://github.com/wayuer19/HyRes_GPU/tree/main/run_OpenMM))  
