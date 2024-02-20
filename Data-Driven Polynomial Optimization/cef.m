function [A,basis] = cef(A,tol)
% The instruction
%
%  [E,BASIS] = CEF(A)
%
% computes a column echelon form E of matrix A
% and returns basis row indices in vector BASIS
%
% The relative threshold for column pivoting can be specified
% as an additional input argument
%
% The reduction is performed by Gaussian elimination with
% column pivoting, based on Matlab's RREF routine

[n,m] = size(A);

% Loop over the entire matrix.
i = 1; j = 1; basis = [];
while (i <= m) & (j <= n)
    % Find value and index of largest element in the remainder of row j
    [p,k] = max(abs(A(j,i:m))); k = k+i-1;
    if (p <= tol)
        % The row is negligible, zero it out
        A(j,i:m) = zeros(1,m-i+1,1);
        j = j + 1;
    else
        % Remember row index
        basis = [basis j];
        % Swap i-th and k-th columns
        A(j:n,[i k]) = A(j:n,[k i]);
        % Find a non-negligible pivot element in the column
        found = 0;
        while ~found
            if abs(A(j,i)) < tol*max(abs(A(:,i)))
                j = j + 1;
                found = (j == n);
            else
                found = 1;
            end;
        end;
        if j <= n,
            % Divide the pivot column by the pivot element
            A(j:n,i) = A(j:n,i)/A(j,i);
            % Subtract multiples of the pivot column from all the other columns
            for k = [1:i-1 i+1:m]
                A(j:n,k) = A(j:n,k) - A(j,k)*A(j:n,i);
            end
            i = i + 1;
            j = j + 1;
        end
    end
end