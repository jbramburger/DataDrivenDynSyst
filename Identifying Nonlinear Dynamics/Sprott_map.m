% -------------------------------------------------------------------------
% Learning Parameter-Dependent Poincare Maps in the Sprott System
%
% We use the third-order chaotic ordinary differential equation:
%
%        x''' + mu*x'' - (x')^2 + x = 0,
%
% where mu is the bifurcation parameter. The user is directed to the paper
% "Simplest dissipative chaotic flow" by J.C. Sprott (Phys. Lett. A, 228,
% pp. 271-274, 1997) for a complete discussion of the system dynamics and 
% its bifurcations. 
%
% TBD...
%
% This script accompanies Section 3.3 of Data-Driven Methods for
% Dynamic Systems. 
%
% For applications of this method to other systems see the repository:
%      https://github.com/jbramburger/Stabilizing_UPOs
%
% Author: Jason J. Bramburger
% -------------------------------------------------------------------------

% Clean workspace
clear all; close all; clc

%% 