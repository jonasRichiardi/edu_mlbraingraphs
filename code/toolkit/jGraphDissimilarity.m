function d=jGraphDissimilarity(g1,g2,varargin)
% JGRAPHDISSIMILARITY computes dissimilarity between two graphs that
%   have fixed vertex ordering (unique node labels) and fixed vertex
%   set cardinality
%
% USAGE:
%   d=jGraphDissimilarity(g1,g2);
%
% INPUT:
%   g1: nV x nV matrix of real numbers - weighted adjacency matrix for graph 1
%   g2: ditto
%   'dissMeasure': string, which dissimilarity measure to use. See code
%       below for choices
%   'PELDcost':  cost for missing edges in the Penalised Edge Label
%       Distance measure. Should be higher than max expected differences
%       (modality dependent).
%
% OUTPUT:
%   d: scalar - a measure of dissimilarity of the graphs
%
% REFERENCE
% J. Richiardi, D. Van De Ville, K. Riesen, H. Bunke, Vector space
%   embedding of undirected graphs with fixed-cardinality vertex sequences
%   for classification, Proc. 2010 Int. Conf. on Pattern Recognition
% 
% VERSION
% v1.0 Dec 2009 Jonas Richiardi
% - initial release, based on discussions with Kaspar Riesen
% v1.1 Aug 2011 JR
% - support for other measures than ICPR 2010 sum-of-edge label diffs
% v1.1.1 July 2013 JR
% - added one-norm for completeness (special case of PELD with no missing
% edges)

[DEBUGMODE,dissMeasure,PELDcost]=process_options(varargin,'DEBUGMODE',false,...
    'dissMeasure','PELD','PELDcost',3);

%% sanity check (can enable here if not running graphPreprocessAdjmat)
% check each adjmat for sanity
% check1=isSaneAdjacencyMatrix(g1);
% check2=isSaneAdjacencyMatrix(g2);
% if ~isempty(check1)
%     disp(check1);
%     error('G1 not an adjacency matrix');
% end
% if ~isempty(check2)
%     disp(check2);
%     error('G2 not an adjacency matrix');
% end
% check graphs have same number of vertices
% if ~all(size(g1)==size(g2))
%     error('Graphs must have same number of vertices');
% end

%% COMPUTE DISTANCE
switch dissMeasure
    case 'PELD'
        %convert upper-triangular adjacency matrix to a vector
        g1_v=jUpperTriMatToVec(g1,1);
        g2_v=jUpperTriMatToVec(g2,1);
        % compute pairwise distance
        d=abs(g1_v-g2_v);
        % add large penalty term (>max possible dist) for missing edges
        maxCostLidx=xor(g1_v~=0,g2_v~=0); % missing edges in ONE of the graphs
        d(maxCostLidx)=PELDcost;
        % total distance is sum of distances
        d=sum(d);
    case 'FrobeniusNorm'
        d=norm(g1-g2,'fro');
    case 'infinityNorm'
        d=norm(g1-g2,inf);
    case 'oneNorm'
        d=norm(g1-g2,1);
    case 'twoNorm'
        d=norm(g1-g2,2);
    case 'kyFanNorm'
        [U,S,V]=svd(g1-g2);
        dS=diag(S);
        d=sum(dS(1:end)); % when using all, is the trace Norm
    case 'SchattenNorm'
        [U,S,V]=svd(g1-g2);
        dS=diag(S);
        p=3;    % if p=2: Frobenius; if p=1; trace norm
        d=(sum(dS.^p)).^(1/p);
    case 'cosNorm'
        g1_v=jUpperTriMatToVec(g1,1);
        g2_v=jUpperTriMatToVec(g2,1);
        d=1-(g1_v/norm(g1_v))'*(g2_v/norm(g2_v));
    otherwise
        error(['Unknown dissimilarity measure ' dissMeasure]);
end
