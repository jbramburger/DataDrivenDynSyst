# Data-Driven Methods for Dynamic Systems

This repository contains the scripts and notebooks that accompany the book Data-Driven Methods for Dynamic Systems. 

The goal of this textbook is to provide an example-driven understanding of how modern computational tools can be applied to interpret dynamic data. In particular, the methods draw inspiration from problems and techniques in the theory of dynamical systems. Sections are organized by starting with a problem, theoretical technique, or concept from dynamical systems and then demonstrating data-driven computational methods based on it. Chapters are organized around a central theme (see below), while sections seek to demonstrate a particular method and the theory it is evoking.

The data that is produced to exemplify each of the methods throughout is generated according to an explicitly stated dynamical system. Examples will focus on some of the most well-studied dynamical systems in the literature, which therefore means we will have access to exact solutions and answers to our questions to gauge the performance of the methods. It should be noted that these example systems are only used to generate data, while methods are agnostic of such an underlying rule unless it is required to implement it. The methods encountered throughout the book are truly data-driven in that they need only be seeded with data from a system, but not information about the system that generates the data. Thus, one may readily apply these methods to dynamic data of ones choice with the examples throughout this book acting as proofs-of-concept. 

## **Packages and Versions**

MATLAB scripts were originally written and run on version R2019_a. Jupyter notebooks use Python 3.8.0 and neural networks are built and trained using TensorFlow 2.0. The reader is referred to the [TensorFlow tutorial](https://github.com/instillai/TensorFlow-Course) to familiarize themself. 

Scripts related to work in Chapter 4 (Data-Driven Polynomial Optimization) require YALMIP and MOSEK to run. Both packages can be download for free at 
- YALMIP: https://yalmip.github.io/download/
- MOSEK: https://www.mosek.com/downloads/

To improve numerical performance, some work related to scripts in Chapter 4 use a Chebyshev function basis instead of monomials. This requires the Chebfun package for MATLAB, which can be freely downloaded at: https://www.chebfun.org/download/

## **Repository Contents**
This repository currently contains four folders, each associated to a chapter of the text. They are organized as follows:

- **Chapter 1**: [**Dynamical Systems: Old and New**](https://github.com/jbramburger/DataDrivenDynSyst/tree/main/Dynamical%20Systems%20Old%20and%20New). This folder contains MATLAB scripts to reproduce the results from Chapter 1 of the textbook. More details can be found in the heading of each script. Organization is as follows:
    - spiral_POD.m applies proper orthogonal decomposition (POD) to a spiral wave solution of a PDE.

- **Chapter 2**: [**Linear Evolution Models**](https://github.com/jbramburger/DataDrivenDynSyst/tree/main/Linear%20Evolution%20Models). This folder contains MATLAB scripts to reproduce the results from Chapter 2 of the textbook. More details can be found in the heading of each script. Organization is as follows:
    - `DMD_Schrondinger.m` applies dynamic mode decomposition (DMD) to the Schrodinger PDE.
    - `windowDMD.m` applies the [windowed DMD method of Dylewsky et al.](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.99.063311) to a multiscale signal to extract the fast and slow timescale dynamics. The MATLAB data `ortho_mat.mat` contains the exact orthogonal matrix used to mix the signal so that one can fully reproduce the figures and results from the textbook.
    - `DelayDMD.m` uses delay coordinates and the Hankel matrix to inflate the dimension of observed data to apply DMD.
    - `EDMD.m` implements [Extended Dynamic Mode Decomposition](https://link.springer.com/article/10.1007/s00332-015-9258-5) to approximate the action of the Koopman operator on the span of observable functions using only data gathered from the system.
    - `Kernel_DMD.m` uses the kernel trick to identify Koopman eigenfunctions from data, based on the work of [Williams et al](https://www.aimsciences.org/article/doi/10.3934/jcd.2015005).  

- **Chapter 3**: [**Identifying Nonlinear Dynamics**](https://github.com/jbramburger/DataDrivenDynSyst/tree/main/Identifying%20Nonlinear%20Dynamics). This folder contains MATLAB scripts to reproduce the results from Chapter 3 of the textbook. More details can be found in the heading of each script. Organization is as follows:
    - `SINDy.m` implements the [sparse identification of nonlinear dynamics method (SINDy)](https://www.pnas.org/doi/10.1073/pnas.1517384113) for nonlinear system identication from data. Also includes implementation the weak formulation due to [McCala and Schaeffer](https://pubmed.ncbi.nlm.nih.gov/28950639/) which is more robust to noisy data.
    - `SINDy_map.m` applies the SINDy method to the discovery of Poincare maps, as outlined by [Bramburger and Kutz](https://www.sciencedirect.com/science/article/pii/S0167278919305470). Particular attention is drawn to the effect of the sparisty parameter and the library.
    - `deJong_control.m` and `Sprott_control.m` implements the [control of chaos method](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.64.1196) on a map and a continuous ODE, respectively. 
    - averaging.m coarse-grains a multiscale signal (using window DMD) and implements the SINDy method to learn the slow-timescale evolution. Based on the work of [Bramburger, Dylewsky, & Kutz](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.102.022204). This script uses the `jupiter_data.mat` and `saturn_data.mat` data files which provide the timeseries data of the position of Jupiter and Saturn in a Sun-Jupiter-Saturn three-body problem.
    - `conserved_quantities.m` implements the method [Kaiser et al.](http://eurika-kaiser.com/downloads/KaKuBr2018cdc.pdf) to learn conserved quantities from data 

- **Chapter 4**: [**Data-Driven Polynomial Optimization**](https://github.com/jbramburger/DataDrivenDynSyst/tree/main/Data-Driven%20Polynomial%20Optimization). This folder contains MATLAB scripts to reproduce the results from Chapter 4 of the textbook. More details can be found in the heading of each script. [YALMIP](https://yalmip.github.io/download/) and [MOSEK](https://www.mosek.com/downloads/) are required to run most scripts in this folder. Organization is as follows:
    - `MG_LyapFn.m` uses polynomial optimization and sum-of-squares relaxations to identify a Lyapunov function for the Moore-Greitzer system.
    - `heteroclinic_UpperBnd.m` and `heteroclinic_LowerBnd.m` identifies an auxialary function to construct barriers in phase space that prove existence and non-existence, respectively, of a heteroclinic orbit in a planar dynamical system.
    - `Disc_LyapFn_Data.m` applies an EDMD-type process to approximate the Koopman operator in order to learn Lyapunov functions from data, as described by [Bramburger and Fantuzzi](https://arxiv.org/abs/2303.01483).
    - `logistic_bounds.m` applies an EDMD-type process to bound expectations long-time averages in a stochastic logistic map.
    - `pend_control.m` discovers a controller from data to stabilize the inverted pendulum on a cart in the upright position. 
    - `invariant_measure.m` identifies extremal invariant measures from data. The method combines the data-driven approximation of the Lie derivative with the method of convex computation of invariant measures due to [Korda et al.](https://arxiv.org/abs/1807.08956)

- **Chapter 5**: [**Learning Dynamics with Neural Networks**](https://github.com/jbramburger/DataDrivenDynSyst/tree/main/Learning%20Dynamics%20with%20Neural%20Networks). This folder contains MATLAB scripts and Jupyter notebooks to reproduce the results from Chapter 5 of the textbook. More details can be found in the heading of each script or notebook. This folder also contains trained neural networks that can be loaded in to reproduce the results from the textbook. Organization is as follows:
    - `neural_network.m` implements a basic neural network and gradient descent training process.
    - `Forecast.ipynb` uses a neural network, implemented and trained using Tensorflow 2.0, to forecast the dynamics of the Henon mapping.
    - `Diffusion_PINN.ipynb` simulations the solutions to the heat equation using a [physics informed neural network (PINN)](https://www.sciencedirect.com/science/article/pii/S0021999118307125). PINN code is repurposed with permission from git user [janblechschmidt](https://github.com/janblechschmidt/PDEsByNNs).
    - `Bistable_PINN.ipynb` uses a PINN to identify the speed of traveling waves in a bistable reaction-diffusion equation.
    - `Diffusion_Discovery.ipynb` employs a PINN to learn the coefficients in a heat equation from data. 

- **Chapter 6**: [**Autoencoder Neural Networks**](https://github.com/jbramburger/DataDrivenDynSyst/tree/main/Autoencoder%20Neural%20Networks). This folder contains MATLAB scripts and Jupyter notebooks to reproduce the results from Chapter 6 of the textbook. More details can be found in the heading of each script or notebook. This folder also contains trained neural networks that can be loaded in to reproduce the results from the textbook. Organization is as follows:
    - `Tent2Logistic.ipynb` and `Tent2Sine.ipynb` both use an autoencoder neural network structure to approximate the conjugacy between the tent map and the logistic and sine maps, respectively.
    - `NormalForm_sn.ipynb` and `NormalForm_pd.ipynb` learn changes of variable to topologically equivalent normal forms in the neighbourhood of codimension 1 bifucations.
    - `GlobalLinearization.ipynb` uses an autoencoder to learn Koopman eignefunctions from data, following [Lusch et al.](https://www.nature.com/articles/s41467-018-07210-0)
    - `ActionAngle.ipynb` learns an invertible change of variable to put the Kepler problem in action-angle coordinates. Much of this code was provided to me by Bethany Lusch and Craig Gin for which I am very thankful. Credit and links to their own GitHub profiles are provided in the corresponding notebook.
    - `Rossler_conj.ipynb` and `Gissinger_conj.ipynb` combine dimensionality reduction and model discovery to identify conjugate maps of their Poincare map dynamics. This work is based on the work of [Champion et al.](https://www.pnas.org/doi/abs/10.1073/pnas.1906995116) and primarily follows [Bramburger et al.](https://www.sciencedirect.com/science/article/pii/S0167278921001652)
