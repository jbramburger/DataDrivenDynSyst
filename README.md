# Data-Driven Methods for Dynamic Systems

MORE INFO COMING SOON - please contact me for details if you are interested.

This repository contains the scripts and notebooks that accompany the book Data-Driven Methods for Dynamic Systems

## **Versions**

MATLAB scripts were originally written and run on version R2019_a. Jupyter notebooks use Python 3.8.0 and neural networks are built and trained using TensorFlow 2.0. The reader is referred to the [TensorFlow tutorial](https://github.com/instillai/TensorFlow-Course) to familiarize themself. 

## **Repository Contents**
This repository currently contains four folders, each associated to a chapter of the text. They are organized as follows:

- **Chapter 1**: [**Dynamical Systems: Old and New**](https://github.com/jbramburger/DataDrivenDynSyst/tree/main/Dynamical%20Systems%20Old%20and%20New). This folder contains MATLAB scripts to reproduce the results from Chapter 1 of the textbook. More details can be found in the heading of each script. Organization is as follows:
    - spiral_POD.m applies proper orthogonal decomposition (POD) to a spiral wave solution of a PDE.

- **Chapter 2**: [**Linear Evolution Models**](https://github.com/jbramburger/DataDrivenDynSyst/tree/main/Linear%20Evolution%20Models). This folder contains MATLAB scripts to reproduce the results from Chapter 2 of the textbook. More details can be found in the heading of each script. Organization is as follows:
    - DMD_Schrondinger.m applies dynamic mode decomposition (DMD) to the Schrodinger PDE.
    - windowDMD.m applies the [windowed DMD method of Dylewsky et al.](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.99.063311) to a multiscale signal to extract the fast and slow timescale dynamics. The MATLAB data ortho_mat.mat contains the exact orthogonal matrix used to mix the signal so that one can fully reproduce the figures and results from the textbook.
    - DelayDMD.m uses delay coordinates and the Hankel matrix to inflate the dimension of observed data to apply DMD.
    - EDMD.m implements [Extended Dynamic Mode Decomposition](https://link.springer.com/article/10.1007/s00332-015-9258-5) to approximate the action of the Koopman operator on the span of observable functions using only data gathered from the system.
    - Kernel_DMD.m uses the kernel trick to identify Koopman eigenfunctions from data, based on the work of [Williams et al](https://www.aimsciences.org/article/doi/10.3934/jcd.2015005).  

- **Chapter 3**: [**Identifying Nonlinear Dynamics**](https://github.com/jbramburger/DataDrivenDynSyst/tree/main/Identifying%20Nonlinear%20Dynamics). This folder contains MATLAB scripts to reproduce the results from Chapter 3 of the textbook. More details can be found in the heading of each script. Organization is as follows:
    - SINDy.m implements the [sparse identification of nonlinear dynamics method (SINDy)](https://www.pnas.org/doi/10.1073/pnas.1517384113) for nonlinear system identication from data. Also includes implementation the weak formulation due to [McCala and Schaeffer](https://pubmed.ncbi.nlm.nih.gov/28950639/) which is more robust to noisy data.
    - SINDy_map.m applies the SINDy method to the discovery of Poincare maps, as outlined by [Bramburger and Kutz](https://www.sciencedirect.com/science/article/pii/S0167278919305470). Particular attention is drawn to the effect of the sparisty parameter and the library.
    - deJong_control.m and Sprott_control.m implements the [control of chaos method](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.64.1196) on a map and a continuous ODE, respectively. This script uses the jupiter_data.mat and saturn_data.mat data files which provide the timeseries data of the position of Jupiter and Saturn in a Sun-Jupiter-Saturn three-body problem. 
    - averaging.m coarse-grains a multiscale signal (using window DMD) and implements the SINDy method to learn the slow-timescale evolution. Based on the work of [Bramburger, Dylewsky, & Kutz](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.102.022204).
    - conserved_quantities.m implements the method [Kaiser et al.](http://eurika-kaiser.com/downloads/KaKuBr2018cdc.pdf) to learn conserved quantities from data 

- **Chapter 4**: [**Polynomial Optimization**]. This folder contains MATLAB scripts to reproduce the results from Chapter 4 of the textbook. More details can be found in the heading of each script. [YALMIP](https://yalmip.github.io/download/) and [MOSEK](https://www.mosek.com/downloads/) are required to run most scripts in this folder. Organization is as follows:
    - TBD

- **Chapter 5**: [**Learning Dynamics with Neural Networks**](https://github.com/jbramburger/DataDrivenDynSyst/tree/main/Learning%20Dynamics%20with%20Neural%20Networks). This folder contains MATLAB scripts and Jupyter notebooks to reproduce the results from Chapter 5 of the textbook. More details can be found in the heading of each script or notebook. This folder also contains trained neural networks that can be loaded in to reproduce the results from the textbook. Organization is as follows:
    - neural_network.m implements a basic neural network and gradient descent training process.
    - Forecast.ipynb uses a neural network, implemented and trained using Tensorflow 2.0, to forecast the dynamics of the Henon mapping.
    - Diffusion_PINN.ipynb simulations the solutions to the heat equation using a [physics informed neural network (PINN)](https://www.sciencedirect.com/science/article/pii/S0021999118307125). PINN code is repurposed with permission from git user [janblechschmidt](https://github.com/janblechschmidt/PDEsByNNs).
    - Bistable_PINN.ipynb uses a PINN to identify the speed of traveling waves in a bistable reaction-diffusion equation.
    - Diffusion_Discovery.ipynb employs a PINN to learn the coefficients in a heat equation from data. 

- **Chapter 6**: [**Autoencoder Neural Networks**](https://github.com/jbramburger/DataDrivenDynSyst/tree/main/Autoencoder%20Neural%20Networks). This folder contains MATLAB scripts and Jupyter notebooks to reproduce the results from Chapter 6 of the textbook. More details can be found in the heading of each script or notebook. This folder also contains trained neural networks that can be loaded in to reproduce the results from the textbook. Organization is as follows:
    - Tent2Logistic.ipynb and Tent2Sine.ipynb both use an autoencoder neural network structure to approximate the conjugacy between the tent map and the logistic and sine maps, respectively.
    - NormalForm_sn.ipynb and NormalForm_pd.ipynb learn changes of variable to topologically equivalent normal forms in the neighbourhood of codimension 1 bifucations.
    - GlobalLinearization.ipynb uses an autoencoder to learn Koopman eignefunctions from data, following [Lusch et al.](https://www.nature.com/articles/s41467-018-07210-0)
    - ActionAngle.ipynb learns an invertible change of variable to put the Kepler problem in action-angle coordinates. Much of this code was provided to me by Bethany Lusch and Craig Gin for which I am very thankful. Credit and links to their own GitHub profiles are provided in the corresponding notebook.
    - Rossler conj.ipynb and Gissinger conj.ipynb combine dimensionality reduction and model discovery to identify conjugate maps of their Poincare map dynamics. This work is based on the work of [Champion et al.](https://www.pnas.org/doi/abs/10.1073/pnas.1906995116) and primarily follows [Bramburger et al.](https://www.sciencedirect.com/science/article/pii/S0167278921001652)
