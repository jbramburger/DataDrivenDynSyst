function K = edmd_with_thresholding(PsiX,PsiY,TOL)

% EDMD method with thresholding to clean small coefficients from the
% Koopman approximation

% Written by Giovanni Fantuzzi

% Useful matrices
A = PsiY * PsiY.';
B = PsiX * PsiY.';
C = PsiX * PsiX.';

% Set up optimization problem
m = size(PsiX, 1);
n = size(PsiY, 1);
K = sdpvar(n, m, 'full');
options = sdpsettings('verbose', 0);
options.mosek.MSK_DPAR_INTPNT_CO_TOL_REL_GAP = 1e-14;
options.mosek.MSK_DPAR_INTPNT_CO_TOL_PFEAS = 1e-14;
options.mosek.MSK_DPAR_INTPNT_CO_TOL_MU_RED = 1e-14;
options.mosek.MSK_DPAR_INTPNT_CO_TOL_DFEAS = 1e-14;
options.mosek.MSK_DPAR_INTPNT_CO_TOL_INFEAS = 1e-14;
OBJ = trace(A) - 2.*trace(K*B) + trace( (K.'*K)*C );

% First solve
optimize([], OBJ, options);
REMOVE = abs( value(K) ) <= TOL;

% Iterate if necessary
CLEAN = any(REMOVE(:));
while CLEAN
    optimize(K(REMOVE)==0, OBJ, options);
    NEW_REMOVE = abs( value(K) ) <= TOL;
    CLEAN = any( NEW_REMOVE(:) - REMOVE(:) );
    REMOVE = NEW_REMOVE;
end
K = value(K);