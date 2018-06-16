function X=graphExtractFeatures(CM3D,mod,classNames,cLabs,C_APPROACH,C_APPROACH_OPT,C_PLOT)
% 
% Extract features from graphs
%
% IN
%   CM3D: struct with fields corresponding to modalities. Each field has
%       two subfields, one per class. The subfields are
%       nRegions x nRegions x nSubjects tensors, where one frontal slab
%       CM3D.mod.class1(:,:,i) is an adjacency matrix for modality mod, 
%       class1, subject i.
%   mod: string - name of modality (field in CM3D)
%   classNames: cell array of strings - names of classes
%   cLabs: 1 x sum(nPatterns) vector of class labels
%   C_APPROACH: string, which approach to use for feature extraction:
%       'direct':           use direct embedding
%       'dissimilarity':    use dissimilarity embedding
%       'vprops':           use vertex properties.
%   C_APPROACH_OPT: string, options for the chosen approach
%   C_PLOT: boolean, whether to plot feature extraction results or not
%
% OUT
%   X:     nFeats x sum(nSubjects) vector of feature vectors
%
% REFERENCES
% Direct embedding: several authors, including
%   - K. Wang et al., "Discriminative analysis of early Alzheimer?s disease
%   based on two intrinsically anticorrelated networks with resting-state
%   fMRI" in Proc. MICCAI 2006, pp. 340?347.
%   - R. C. Craddock et al., "Disease state prediction from resting
%   state functional connectivity". Magn. Reson. Med. 62(6), 2006
%   - J. Richiardi et al.,"Vector space embedding of undirected graphs
%   with fixed-cardinality vertex sequences for classification" in Proc
%   ICPR 2010, pp. 902-905.
% Dissimilarity embedding:
%   - our ICPR paper (same as above)
% Vertex Properties (in a machine learning context): several authors, including
%   -  G. A. Cecchi et al. "Discriminative network models of schizophrenia",
%   in Proc. NIPS 2009, pp. 252?260.
%   - J. Richiardi et al., "Classifying connectivity graphs using graph
%   and vertex attributes", Proc. PRNI 2011
%
% VERSION
% v1.0 July 2013 Jonas Richiardi for Grenoble summer school

% sanity check
if ~isstruct(CM3D)
    error('CM3D must be a struct');
end
if ~isfield(CM3D,mod)
    error(['CM3D has no field named ' mod]);
end

% re-compute useful stuff
nSubjects=hist(cLabs,unique(cLabs));
nClasses=numel(nSubjects);
nVertices=size(CM3D.(mod).(classNames{1}),1);

% do feature extraction
switch C_APPROACH
    case 'direct' % direct
        vsubsetIdx=1:nVertices;         % project to a subset of vertices? 
        mynPairs=nchoosek(size(CM3D.(mod).(classNames{1})(vsubsetIdx,vsubsetIdx,1),1),2);        
        X=zeros(mynPairs,sum(nSubjects),'single');
        for c=1:nClasses
            stOffset=sum(nSubjects(1:c-1));
            for s=1:nSubjects(c)
                embedThis=jUpperTriMatToVec(CM3D.(mod).(classNames{c})(vsubsetIdx,vsubsetIdx,s),1);
                X(:,stOffset+s)=embedThis;
            end
        end
    case 'dissimilarity' % dissimilarity        
        CM3Dfull=cat(3,CM3D.(mod).(classNames{1}),CM3D.(mod).(classNames{2}));
        X=jGraphDissimilarityEmbedding(CM3Dfull(:,:,:),...
            CM3Dfull(:,:,:),'dissMeasure',C_APPROACH_OPT);
    case 'vprops' % vprops
        nProps=4; % how many properties we use
        X=zeros(nVertices*nProps,sum(nSubjects),'single');
        for c=1:nClasses
            fprintf('Extracting vertex properties for class %d\n',c);
            stOffset=sum(nSubjects(1:c-1));
            tic
            % extract all properties in turn for all vertices and subjects
            parfor s=1:nSubjects(c)
                fprintf('%d ',s);
                thisA=abs(CM3D.(mod).(classNames{c})(:,:,s));
                v_s=strengths_und(thisA);           % vertex strength - could do this without BCT...
                v_CC=clustering_coef_wu(thisA);     % vertex clustering coefficient
                v_El=efficiency_wei(thisA,1);       % local efficiency
                v_bw=betweenness_wei(thisA);        % betweenness centrality
                embedThis=[v_s(:);v_CC(:);v_El(:);v_bw(:)];     % concatenate
                X(:,stOffset+s)=embedThis;
            end
            disp('');
            toc
        end
    otherwise
        error(['Unknown feature extraction method ' C_APPROACH])
end
  
% plot if needed
if C_PLOT
    maxVal=max(X(:));
    labelling=(cLabs/max(cLabs))*maxVal;
    labelling=repmat(labelling,ceil(size(X,1)/50),1);
    figure; imagesc([X;labelling]); colorbar;
    title([classNames{1} '(first ' num2str(nSubjects(1)) ')' ...
        'vs. ' classNames{2} '(last ' num2str(nSubjects(2)) ')'])
    xlabel('examples'); ylabel('features');

end