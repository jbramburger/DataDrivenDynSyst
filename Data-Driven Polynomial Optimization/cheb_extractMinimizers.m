function x = cheb_extractMinimizers(momentMatrices, gramMonomials)

% Attempt to extract minimizers from a set of moment matrices. Inputs are:
% - momentMatrices: cell array of moment matrices
% -  gramMonomials: indices of the chebyshev polynomials used in the
%                   construction of the moment matrix

tol = 1e-12;
droptol = 1e-1;
if iscell(momentMatrices)
    for k = 1:length(momentMatrices)
        % Rank via SVD decomposition
        cleantol = tol*max(max(abs(momentMatrices{k})));
        momentMatrices{k} = clean(momentMatrices{k},cleantol);
        [U,S] = svd(momentMatrices{k});
        [S,pos] = sort(diag(S),'descend');
        U = U(:,pos);
        drop = S(2:end)./( eps + S(1:end-1) );
        drop = find(drop<droptol,1,'first');
        if ~isempty(drop)
            rankM = drop;
        else
            rankM = length(S);
        end
        U = U(:,1:drop)*diag(sqrt(S(1:drop)));
        
        % Get column echelon form, like gloptipoly
        [U,basis] = cef(U,1e-6);
        
        % Decompose the powers in v and find monomials of degree 1
        beta = gramMonomials{k}(basis);
        nmons = length(gramMonomials{k});
        N = zeros(rankM^2,1);
        I = []; J = [];
        for j = 1:rankM
            pow = [beta(j)+1; abs(beta(j)-1)];
            [ind,pos] = ismember(pow(:),gramMonomials{k}(:));
            I = [I; j; j];
            J = [J; pos(:)];
        end
        A = sparse(I,J,0.5*ones(size(I)),rankM,nmons);
        N = A*U;
        
        % Order shur decomposition and extract minimizers
        [Q,T] = schur(N);
        
        if isempty(T)
            % complex eigenvalues = failure
            x = [];
            
        else
            % Retrieve optimal vectors
            % It is assumed than there is no multiple root
            x{k} = zeros(1,rankM);
            for i = 1:rankM
                    x{k}(i) = Q(:,i)'*N*Q(:,i);
            end
        end
    end
    
else
    % Make a cell and call again
    x = cheb_extractMinimizers({momentMatrices},{gramMonomials});
    x = x{1};
end







