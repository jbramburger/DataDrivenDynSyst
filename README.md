# Data-Driven Methods for Dynamic Systems

MORE INFO COMING SOON - please contact me for details if you are interested.

This repository contains the scripts and notebooks that accompany the book Data-Driven Methods for Dynamic Systems

Neural networks use TensorFlow 2. A tutorial is available at: https://github.com/instillai/TensorFlow-Course 

## **Repository Contents**
This repository currently contains four folders, each associated to a chapter of the text. They are organized as follows:

- Chapter 2: [**Linear Evolution Models**](https://github.com/jbramburger/DataDrivenDynSyst/tree/main/Linear%20Evolution%20Models): This folder contains MATLAB scripts to reproduce the results from Chapter 2 of the textbook. More details can be found in the heading of each script. Organization is as follows:
    - DMD_Schrondinger.m applies dynamic mode decomposition (DMD) to the Schrodinger PDE.
    - windowDMD.m applies the [windowed DMD method of Dylewsky et al.](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.99.063311) to a multiscale signal to extract the fast and slow timescale dynamics.
    - DelayDMD.m uses delay coordinates and the Hankel matrix to inflate the dimension of observed data to apply DMD.
    - EDMD.m implements [Extended Dynamic Mode Decomposition](https://link.springer.com/article/10.1007/s00332-015-9258-5) to approximate the action of the Koopman operator on the span of observable functions using only data gathered from the system.
    - Kernel_DMD.m uses the kernel trick to identify Koopman eigenfunctions from data, based on the work of [Williams et al](https://www.aimsciences.org/article/doi/10.3934/jcd.2015005)  
