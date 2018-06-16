function vOut=jGraphDissimilarityEmbedding(G,P,varargin)
% JGRAPHDISSIMILARITYEMBEDDING compute dissimilarity embedding of 
% a graph or graphs with respec to a set of prototypes
%
%
% INPUT:
%   G: nVertices x nVertices xnGraphs matrix of floats - weighted adjacency matrix
%       of the graphs to embed
%   P: nVertices x nVertices x nPrototypes 3D matrix of floats -
%       adjacency matrices of prototypes
%
% OUTPUT:
%   vOut: nPrototypes x 1 vector of floats - dissimilarity embedding
%       of G
%
% NOTES:
%   This code works for graphs with fixed-ordering vertex sequences, that
%   is graphs whose vertex set always has the same cardinality, and posess
%   the unique node label property (fixed vertex ordering)
%
% VERSION:
% v1.0 Dec 2009 Jonas Richiardi
%   - initial release
% 1.1 Aug 2011 JR
%   - option to choose dissimilarity measure

[DEBUGMODE,dissMeasure]=process_options(varargin,'DEBUGMODE',false,...
    'dissMeasure','ICPR2010');
nG=size(G,3);
nP=size(P,3);

%% sanity check (can enable here if not running graphPreprocessAdjmat)
% check graphs to embed
% for gidx=1:nG
%     checkG=isSaneAdjacencyMatrix(G(:,:,gidx));
%     if ~isempty(checkG)
%         disp(checkG);
%         error(['Graph ' num2str(gidx) ' is not an adjacency matrix']);
%     end
% end
% % check prototypes
% for pidx=1:nP
%     checkP=isSaneAdjacencyMatrix(P(:,:,pidx));
%     if ~isempty(checkP)
%         disp(checkP);
%         error(['Prototype ' num2str(pidx) ' is not an adjacency matrix']);
%     end
% end

% check graphs have same number of vertices as prototypes
if ~all(size(G(:,:,1))==size(P(:,:,1)))
    error('Graph to embed must have same number of vertices as prototypes');
end

%% do embedding
vOut=zeros(nP,nG);
parfor gidx=1:nG
    for pidx=1:nP
        vOut(pidx,gidx)=jGraphDissimilarity(G(:,:,gidx),P(:,:,pidx),...
            'dissMeasure',dissMeasure);
    end
end
