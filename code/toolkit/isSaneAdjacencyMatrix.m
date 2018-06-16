function status=isSaneAdjacencyMatrix(A,varargin)
% ISSANEADJACENCYMATRIX check whether the adjacency matrix is correctly
% defined
%
% USAGE:
%   status=isSaneAdjacencyMatrix(A,varargin)
% v1.0 Oct 2009 Jonas Richiardi
% - initial release
% v1.1 April 2013 JR
% - faster code
% - could add special code for sparse mats (nnz(A-A.')>0)

[undirected]=process_options(varargin,'undirected',true);

status=[];

% check there are no nans
if any(isnan(A(:)))
    status=[status 'Adjacency matrix must not contain NANs. '];
end

% check that adjmat is square
if (size(A,1)~=size(A,2))
    status=[status 'Adjacency matrix must be square. '];
else
    % check that the adjmat is symmetric
    %if ~(all(all(A==A')))
    if ~isequal(A,A.')
        status=[status 'Adjacency matrix must be SYMMETRIC'];
        problemLocation=mat2str(find(any(A-A.')));
        status=[status ' (problem may be at vertices ' problemLocation ...
            ', although maybe NaNs are to blame). '];
    else
        % check that there are no loops (2*nLoops = tr(g))
        if (undirected==true) && (trace(A)~=0)
            status=[status 'Undirected graph must be ACYCLIC ' ... 
                '- check adjacency matrix for loops. '];
        end
    end
end
