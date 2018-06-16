function graphClassification_cheap
%
% *INTRO*
% This is a simple educational script for
% classifying brain graphs*. It comprises graph preprocessing, feature
% extraction, feature selection, cross-validation, training, testing,
% classifier selection, and information mapping.
% 
% It will not provide state of the art performance, but can be used as an
% exploratory basis for your own work. Apart from modern embeddings, kernels
% and graph CNNs, it lacks nested cross-validation for classifier parameter
% setting and projection dimension in feature selection. This is not done
% here to save time given the short time frame available for teaching this...
%
% move between script cells by using the cell menu in the Matlab editor
%
% Feedback, bug reports, PRs are welcome.
% 
% *REFERENCES*
% * For background, theory, and context, please see Richiardi et al.,
% Machine Learning with Brain Graphs: Predictive Modeling Approaches for
% Functional Imaging in Systems Neuroscience, IEEE Signal Processing
% Magazine, May 2013 (available at http://www.richiardi.net/papers/2013-Richiardi-IEEESigProcMag.pdf)
%
% *ACKNOWLEDGEMENTS*
% S. Achard (INPG), E. Bullmore (Cambridge), H. Bunke (U. Bern),
% J.-F Coeurjolly (INPG), N. Leonardi (EPFL), D. Van De Ville (EPFL),
% K. Riesen (U. Bern)
% Funding: Marie Curie IOF #299500 "Modelling and inference in brain
% networks for diagnosis"
% 
% *VERSION HISTORY*
% v1.0 July 2013 Jonas Richiardi
%   First presented at the Summer school on graphical models for the
%   characterisation of information flow in complex networks: application
%   in neuroimaging
% v1.0.1 Oct 2017
%   Slight modifications for the Verona summer school on brain connectomics
% v1.0.2 June 2018
%   Further small tweaks for the advanced MRI course in the lemanic
%   neuroscience doctoral school
% v1.0.3 July 2018
%   More cleanup for HBM2018 Educational course "pattern recognition for
%   neuroimaging"

%% ENVIRONMENT CHECKING

% add paths to tools we need
% This also has the nice side effect of shadowing the bioinfo toolbox's svmtrain
addpath(genpath('toolkit'));

% check paths - compiled and toolbox classifiers
disp('Checking your environment');
if strcmp(computer('arch'),'glnxa64')
    warning('RF and SVM classifiers have been compiled for Linux but not tested seriously.');
end
% SVM
if exist('svmtrain','file')~=3
    warning('libSVM training function binary not found on path. SVM training will not work. ');
end
if exist('svmpredict','file')~=3
    warning('libSVM prediction function binary not found on path. SVM prediction will not work. ');
end
% RF
if exist('mexClassRF_train','file')~=3
    warning('RandomForest-matlab training function binary not found on path. RF training will not work. ');
end
if exist('mexClassRF_predict','file')~=3
    warning('RandomForest-matlab testing function binary not found on path. RF prediction will not work. ');
end
% NB
if exist('NaiveBayes','class')~=8
    tbcheck=ver('stats');
    if isempty(tbcheck)
        warning('Looks like you do not have the statistics toolbox. Na?ve Bayes will not work.');
    else
        warning('Looks like you have the statistics toolbox but Na?ve Bayes is not found on path. Na?ve Bayes will not work.');
    end
end

% now check for parallel computing - skip this in 'cheap' version
% if exist('matlabpool','file')==2
%     if matlabpool('size')==0
%         disp('Attempting to start parallel (multicore) computing pool');
%         matlabpool('open');
%         if matlabpool('size')==0
%             error('Cannot start matlab pool. This will require some manual editing. Please ask for help.');
%         else
%             disp('It seems we successfully started a matlab pool.');
%         end
%     else
%         disp('Parallel (multicore) pool seems to be already running, good.');
%     end
% else
%     error(['You do not seem to have access to parallel computing. '...
%         'This requires editing a few files manually. Please ask for help.']);
%     % simply replace the "parfor" loop by a "for" loop in these three
%     % files:
%     % graphExtractFeatures.m, graphSelectFeatures.m,
%     % jGraphDissimilarityEmbedding.m
% end

disp('If you have no warnings above the environment looks OK!');
disp('Warnings about pre-existing matlabpool jobs are not an issue.');



%% DATA LOADING
% 
% Choose the dataset. Here only 'UCLA Autism' is provided (both fMRI and DTI data, diagnosis 
% prediction), so set |C_DATASET| to |UCLA_Autism|.
% 
% Also, choose the modality with |C_MOD| where 1 is fMRI, and 2 is DTI.
%
% UCLA Autism dataset: stored in a struct called |CM3D|. fMRI data
% is in the struct |CM3D.fMRI|, and DTI data is in the struct |CM3D.DTI|.
% Within each, the graphs for healthy controls (typically developing)
% subjects is in the field |TD|, and graphs for autistic subjects are in
% the field |ASD|. Each field contains a 3D matrix of dimensions
% nRegions x nRegions x nGraphs.
% 
% Other datasets are stored using a similar principle.
%
% Example 1: to extract the second healthy control's fMRI graph, do this:
%   |myGraph=squeeze(CM3D.fMRI.TD(:,:,2));|
%
% Example 2: to extract and show the 5th ASD subjects' DTI graph, do this:
%   |myGraph=squeeze(CM3D.DTI.ASD(:,:,5));|
%   |figure('name','DTI/ASD subject 5'); imagesc(myGraph); axis square;| 
%   |colorbar; xlabel('vertex index'); ylabel('vertex index');|
%
% Details on experimental paradigm, acquisition, preprocessing, and
% papers where these data appeared are in |data/README.txt|
%
% Also, we load a codebook each time. this is a useful structure called
% |CB|. this has fields
%
% * 'id'        The ID number of the atlas region (intensity in Atlas)
% * 'num'       The number of atlas regions
% * 'name'      Full name of the atlas region
% * 'sname'     Shortened name of the atlas region
% * 'rname'     The name
% * 'center'    MNI coordinates of the center of this atlas region in the
%                   *structural* image
%
% Have a look to get a feel for the data! For example scatter-plot vertex
% centers (MNI) to get a feel for spatial sampling in various datasets.

C_DATASET='UCLA_Autism';
C_MOD=1;           % modality choice (for UCLA_Autism only). 1: fMRI, 2: DTI
C_BALANCECLASSES=true;

switch C_DATASET
    case 'UCLA_Autism'     
        mods={'fMRI','DTI'};                        % modality names
        classNames={'TD','ASD'};                            % unique class labels
        load(fullfile('..','data','UCLA_Autism','CM3D.mat'));
         load(fullfile('..','data','UCLA_Autism','codebook.mat'));
    otherwise
        error('Unknown dataset');
end

nClasses=numel(classNames);                           % number of classes
nGraphs=zeros(1,nClasses);                            % number of patterns per class
nVertices=size(CM3D.(mods{C_MOD}).(classNames{1}),1); % number of regions
cLabs=[];                                             % class labels per sample

% compute number of graphs per class and assign class labels
for c=1:nClasses
    nGraphs(c)=size(CM3D.(mods{C_MOD}).(classNames{c}),3);
    cLabs=[cLabs ones(1,nGraphs(c))*(c-1)];
end

if C_BALANCECLASSES==true
    % check if we have differing number of samples per class
    if (numel(unique(nGraphs))>1)
        cLabs=[];
        minNsamples=min(nGraphs);
        for c=1:nClasses
            cLabs=[cLabs ones(1,minNsamples)*(c-1)];
            if nGraphs(c)>minNsamples
                disp(['Modality ' mods{C_MOD} ': Subsampling class ' classNames{c} ' from start']);
                CM3D.(mods{C_MOD}).(classNames{c})=CM3D.(mods{C_MOD}).(classNames{c})(:,:,1:minNsamples);
            end
        end
        nGraphs=repmat(minNsamples,1,nClasses); % update
    end
end

%% GRAPH PREPROCESSING
% First, we need to make sure no self-loops are present. Then, we might
% want to threshold away connections that are 'small' or 'non-significant'
% in some sense. If we want to use graph properties, we need to make sure
% that all graphs have the same edge density (a sufficient condition if the
% fixed-cardinaly vertex sequence property is fulfilled).
%
% Note that it would be a good idea to remove some vertices depending on
% signal quality or signal drop-out. This depends on the MRI machine,
% sequence, number of channels in the head coil, atlas, and segmentation 
% algorithm. Typically areas close to air cavities may have dropout.
% Small atlas regions (e.g. globus pallidus in the AAL) generally have
% very poor signal, which could be a combination of spatial location
% (ventral) and having few voxels.

C_THRESHMETHOD='edgeDensity';   % '' or 'edgeDensity'
C_THRESHVAL=0.9;                % scalar between 0 and 1 or empty matrix ([])
                                % take modality into account when setting
                                % this! fMRI + Pearson pairwise yields a
                                % complete graph, while DTI + deterministic
                                % tracto will most probably not (typically
                                % around 5-10% unthresholded edge density)

CM3D = graphPreprocessAdjmats(CM3D,mods{C_MOD},classNames,...
    C_THRESHMETHOD,C_THRESHVAL);


%% ADJACENCY MATRIX VISUALISATION
% "There is no excuse for not looking at the data!"
% Here we'll look at a few matrices only, but do modify the code to look 
% at all of them...

figure; 
for c=1:nClasses
    rnd_graphs=floor(rand(1,2)*nGraphs(c))+1;            % select 2 graphs randomly
    for g=1:2
        subplot(2,2,g+c*(c-1));
        imagesc(squeeze(CM3D.(mods{C_MOD}).(classNames{c})(:,:,rnd_graphs(g))));
        axis square; colorbar; 
        title([classNames{c} ' ' num2str(rnd_graphs(g))]);
    end
end

%% FEATURE EXTRACTION / GRAPH REPRESENTATION
% We can now choose a representation and classification approach for our
% graph data and extract features. Here we'll explore *direct embedding*
% of edge labels, *dissimilarity embedding*, and *graph/vertex properties*.
%
% Each subject's features are extracted independently of other
% subjects and class labels, so we don't need to include this in a
% cross-validation setting.


C_APPROACH='direct';            % 'direct': direct embedding
                                % 'dissimilarity': dissimilarity embedding
                                % 'vprops': vertex properties
C_APPROACH_OPT='kyFanNorm';     % option for graph representation approach
                                % for dissimilarity embedding, this is the
                                % dissimilarity function. See
                                % jGraphDissimilarity.m for details and
                                % choices
C_PLOTFX=true;                  % whether to plot extraced features or not

X=graphExtractFeatures_cheap(CM3D,mods{C_MOD},classNames,cLabs,...
    C_APPROACH,C_APPROACH_OPT,C_PLOTFX);
fprintf('Extracted %d features and %d graphs\n',size(X,1),size(X,2));
D=size(X,1);                    % number of feature (pre-feature selection)


%% CLASSIFICATION AND CROSS-VALIDATION FOLDS SETUP
% Now we have a vector representation of graphs, we can choose a
% classifier and its options, how many features we want to select, and 
% 


C_D=2000; % nchoosek(nVertices,2);          % how many features to select (this may be modified by the CV loop if too high)

C_CHECKWITHINCLASSVAR=true;         % check within-class variance for 0 and don't select features that are pathological 
                                    % this is useful for generative models
                                    % that fit class-conditional densities
                                    % (e.g. the NB model). Otherwise skip
                                    % for speed

C_PLOTFS=false;                     % whether to display FS plots or not

C_MYCLASSIFIER='RF';               % 'NB': Naive Bayes
                                    % 'SVM': Support Vector Machine
                                    % 'RF': Random Forest

% classifier parameters - NB
C_MYCPARAMS.NB.Distribution='normal'; % which type of density estimate to use
                                    %   'normal', 'kernel'
                                    
% classifier parmeters - SVM
C_MYCPARAMS.SVM.kernel=0;         % 0: linear (u'*v). 1 polynomial, 2. RBF
C_MYCPARAMS.SVM.degree=0;         % degree in polynomial kernel function

                                
% classifier parameters - RF
C_MYCPARAMS.RF.nTrees=501;        % number of trees in a RF
C_MYCPARAMS.RF.nFeats=round((sqrt(C_D)));         % number of features per tree in a RF
C_MYCPARAMS.RF.replace=true;      % RF: sample with replacement
C_MYCPARAMS.RF.nodesize=1;        % RF: minimum number of cases per leaf


%nFolds=sum(nGraphs);                % leave-one-subject-out cross validation
nFolds=10;                           % k-fold cross-val
foldLabels=randsample([repmat(1:nFolds,1,floor(sum(nGraphs)/nFolds)) ...
    1:mod(sum(nGraphs),nFolds)],sum(nGraphs));
%figure; hist(foldLabels,unique(foldLabels)); axis tight; xlabel('fold');
%ylabel('number of left out samples per fold');


%% DO CV folds
fprintf('Starting %d CV folds...\n',nFolds);
allPredictions=[];                  % this wil store classifier predictions
allTruths=[];                       % this will store ground truths
allBestFeatsIdx=cell(nFolds,1);     % cell array to save the selected features at each fold
allCls=cell(nFolds,1);              % cell array to save the trained classifiers at each fold
allTRidx=cell(nFolds,1);            % save training sample indices
allcLabsTR=cell(nFolds,1);          % save training labels

tic
for f=1:nFolds
    fprintf('%d ',f);
    %clear X_TR X_TE cLabs_TR cLabs_TE;
    X_TR=[]; X_TE=[]; cLabs_TR=[]; cLabs_TE=[];

    %%% Separate training and testing set
    % This is essential for proper error estimation.
    TRsIdx=foldLabels~=f;   % take al data who is not tagged to be left out this fold
    TEsIdx=~TRsIdx;
    TRdataIdx=TRsIdx;   % for expansion
    allTRidx{f}=TRdataIdx;
    TEdataIdx=TEsIdx;   % for expansion
    if strcmp(C_APPROACH,'dissimilarity')
        X_TR=X(TRdataIdx,TRdataIdx);
        X_TE=X(TRdataIdx,TEdataIdx);
    else
        X_TR=X(:,TRdataIdx); % cut out training data
        X_TE=X(:,TEdataIdx); % cut out testing data
    end
    cLabs_TR=cLabs(TRdataIdx);
    allcLabsTR{f}=cLabs_TR;
    cLabs_TE=cLabs(TEdataIdx);
    
    %%% FEATURE SELECTION
    % Feature selection is trained on the training set only, by seeing
    % class labels, then applied to both the train and test sets.
    
    if C_D>size(X_TR,1)
        warning(['C_D is too large at ' num2str(C_D) ...
            ', adjusting and using all features (' num2str(size(X_TR,1)) ')']);
        C_D=size(X_TR,1);
        bestFeaturesIdx=[1:size(X_TR,1)]';
    elseif C_D==size(X_TR,1) % no feat sel
        bestFeaturesIdx=[1:size(X_TR,1)]';
    else
        bestFeaturesIdx=graphSelectFeatures(X_TR,cLabs_TR,C_D,...
            C_CHECKWITHINCLASSVAR,C_PLOTFS);

        % To see the dramatic impact of cheating, uncomment next line - this allows the
        % feature selection procedure to see the whole dataset and labels.
        % bestFeaturesIdx=graphSelectFeatures(X,cLabs,C_D,C_CHECKWITHINCLASSVAR,C_PLOTFS);
    end
    allBestFeatsIdx{f}=bestFeaturesIdx;
    X_TR_FS=X_TR(bestFeaturesIdx,:);
    X_TE_FS=X_TE(bestFeaturesIdx,:);
    
    %%% normalise if needed
    % SVMs are sensitive to input scales (as opposed to NB or RF), so if
    % we're using a representation where scales can differ importantly
    % (e.g. in vertex properties embedding), performance can degrade.
    % to solve this, we will simply rescale each feature to unit std.
    % we could also z-norm, or use other strategies
    if strcmp(C_MYCLASSIFIER,'SVM')
        % estimate std on training set (don't look at test set!)
        sigma_TR_FS=std(X_TR_FS,[],2);
        % normalise training set
        X_TR_FS=X_TR_FS./repmat(sigma_TR_FS,1,size(X_TR_FS,2));
        % normalise test set
        X_TE_FS=X_TE_FS./repmat(sigma_TR_FS,1,size(X_TE_FS,2));
    end
    
    %%% TRAINING
    % At this stage we now have feature vectors, potentially reduced
    % dimension through feature selection. We're ready to train the
    % classifier on the training set.
    %
    % For Naive Bayes, we use the statistics toolbox implementation by the 
    % Mathworks.
    %
    % For Support Vector Machines, we use the libSVM implementation by 
    % Chih-Chung Chang and Chih-Jen Lin.
    %
    % For Random Forests, we use the randomforest-matlab implementation by
    % Abhishek Jaiantilal, based on work by Leo Breiman, Adele Cutler,
    % Andy Liaw and Matthew Wiener.
    %clear myCls;
    myCls=[];
    switch C_MYCLASSIFIER
        case 'NB'
            myCls = NaiveBayes.fit(X_TR_FS', cLabs_TR);
        case 'SVM'
            SVM_T=['-t ' num2str(C_MYCPARAMS.SVM.kernel) ' '];
            SVM_D=['-d ' num2str(C_MYCPARAMS.SVM.degree) ' '];            
            myCls=svmtrain(double(cLabs_TR'),double(X_TR_FS'),[SVM_T SVM_D ' -h 1 -c 1']);
        case 'DT'
            myCls = classregtree(X_TR_FS',cLabs_TR,'method','classification');
        case 'RF'
            myCls = classRF_train(double(X_TR_FS'),double(cLabs_TR),...
                C_MYCPARAMS.RF.nTrees,C_MYCPARAMS.RF.nFeats,...
                C_MYCPARAMS.RF);
        otherwise
            error('Unknown classifier');
    end
    allCls{f}=myCls;
    
    %%% TESTING
    % Now the classifier has been trained. We need to see how well it does
    % against unseen data - the held-out testing set.
    
    %clear predLabels
    predLabels=[];
    switch C_MYCLASSIFIER
        case 'NB'
            predLabels = predict(myCls,X_TE_FS');
        case 'SVM'
            [predLabels,~,predVals]=svmpredict(double(cLabs_TE'),double(X_TE_FS'),myCls);
        case 'DT'
            [predLabels,~] = myCls.eval(X_TE_FS');
            predLabels=cell2num(predLabels);
        case 'RF'
            predLabels = classRF_predict(double(X_TE_FS'),myCls);
        otherwise
            error('Unknown classifier');
    end
    % save predictions and truths for later analysis
    allPredictions=[allPredictions; predLabels(:)];
    allTruths=[allTruths;cLabs_TE(:);];
    
end % CV folds
disp('Done with Cross-validation.');
toc

% compute and print overall and class accuracies for convenience
myStats=jEvaluateClassifierNClassesDiscreteLabels(allTruths+1,allPredictions+1,'nClasses',2);
disp(' acc   c0   c1');
disp(mat2str([myStats.meanClassAcc myStats.classAccs],2));
myStats.confMat

%% VISUALISATION OF FEATURE SPACE AND DECISION BOUNDARY
% We can have a look at the training data set in feature space by using
% the two features with the highest point-biserial correlation with the
% class.
% 
%
% Please set |C_D=2| for proper operation. The decision boundary will be 
% meaningless if this setting produces chance results in CV.
% 
figure; hold on; 
plot(X_TR_FS(1,cLabs_TR==0),X_TR_FS(2,cLabs_TR==0),'rx','MarkerSize',10);
plot(X_TR_FS(1,cLabs_TR==1),X_TR_FS(2,cLabs_TR==1),'go','MarkerSize',10);
xlabel(['Feature 1: ' num2str(bestFeaturesIdx(1))]);
ylabel(['Feature 2: ' num2str(bestFeaturesIdx(2))]);
legend(['Class 0: ' classNames{1}],['Class 1: ' classNames{2}]); axis equal;
if C_D==2
    myMap=jPlotDecisionBoundary(myCls,C_MYCLASSIFIER,X_TR_FS);
else
    warning('Not plotting decision boundary, feature space dimensionality > 2');
end
%xlim([-10 10]); ylim([-10 10]);

%% INFORMATION MAPPING
% It is interesting to have a look at which features drive the
% classification. Extraction procedure vary depending on classifier, but
% ultimately results in a map indicating the relative importance of
% features. What the features are depends on the embedding type.
% 
% Direct embedding yields feature importances in edge space.
% Vertex property embedding yields feature importances in vertex space.
% Dissimilarity embedding yields feature importance in graph space.
% 

% gather stats on weight map across folds
myInfoMap=zeros(D,1);
disp('Computing feature importances');
for f=1:nFolds
    switch C_MYCLASSIFIER
        case 'SVM'
            % first check if model is linear. If not, abandon ship.
            if C_MYCPARAMS.SVM.kernel==2
                error('Code not ready to map RBF kernels');
            elseif C_MYCPARAMS.SVM.kernel==0 || (C_MYCPARAMS.SVM.kernel==1 && C_MYCPARAMS.SVM.degree==1)
                % compute primal weight vector. class0=+1, class1=-1, so no need
                % to flip the weight vector.
                thisW=allCls{f}.SVs'*allCls{f}.sv_coef;
                myInfoMap(allBestFeatsIdx{f})=myInfoMap(allBestFeatsIdx{f})+thisW.^2;
            else
                error('Can only map linear kernels');
            end
        case 'RF'
            myInfoMap(allBestFeatsIdx{f})=myInfoMap(allBestFeatsIdx{f})+allCls{f}.importance;
        case 'NB'
            theseImps=jNBFeatureImportance(X(allBestFeatsIdx{f},allTRidx{f}),allcLabsTR{f}+1);
            myInfoMap(allBestFeatsIdx{f})=myInfoMap(allBestFeatsIdx{f})+theseImps(:);
        otherwise
            error('Code not ready for mapping other classifiers');
    end
end

% compute summary maps. Which one to use depends on intent. The normalised
% info map loses the original meaning of feature importance (weight in SVM,
% Gini importance in RF, odds ratio in NB) but still displays relative
% importance and avoids scaling headaches for visualisation...
myInfoMap=myInfoMap/nFolds;             % average importance for RF, mean square for SVM
myInfoMapN=myInfoMap/sum(myInfoMap);    % sum to 1

% output the info map in a sensible way
switch C_APPROACH
    case 'direct'
        % info map is represented as edges of a discriminative graph.
        % Note we could also aggregate all edge importances to their
        % corresponding region and get an approximation of the
        % discriminative power of brain regions.
        A_map=jVecToSymmetricMat(myInfoMapN,size(CM3D.(mods{C_MOD}).(classNames{1}),1));
        % show discriminative graph
        figure; imagesc(A_map); colorbar;
        title(['discriminative graph for ' C_DATASET ' with ' C_MYCLASSIFIER],...
            'Interpreter','none');
        jAdjMatToNetworkFile(A_map,CB.sname,'infoMap.txt','formatOut','cytoScape');
    case 'dissimilarity'
        % info map shows how influential each graph is in the dissimilarity
        % representation
        [mySortMap,mySortMap_idx]=sort(myInfoMapN,'descend');
        figure; plot(mySortMap,'LineWidth',3,'Color','k'); xlabel('graph index');
        ylabel([C_MYCLASSIFIER ' importance']); hold on;
        c=cell(nClasses,1); myLegend=cell(nClasses+1,1);
        myLegend{1}='importance';
        cMarkers={'go','rx'};
        for c=1:nClasses
            cPos{c}=find(cLabs(mySortMap_idx)==c-1);
            scatter(cPos{c},mySortMap(cPos{c}),cMarkers{c},'SizeData',150,'LineWidth',3);
            myLegend{c+1}=classNames{c};
        end
        legend(myLegend); axis tight;
    case 'vprops'
        % info map shows how discriminative a region is with respect to a
        % particular vertex attribute
        [sXbase,sYbase,sZbase] = sphere(30); % coord for sphere -> draw sphere with surf (filled) or mesh (grid)
        myCData=ones(size(sXbase))*1;   % choose color scheme
        sphereScale=300;            % how to linearly scale the spheres (depends on coordinate system, classifier...)
        showTopNlabels=10;          % how many vertex labels to show
        nProps=3;                   % if changing here must also change graphExtractFeatures
        pNames={'strength','clustering coefficient','betweenness centrality'};
        load(fullfile('.','toolkit','convexHullOutlinePlanes.mat'));
        for p=1:nProps
            stOffset=nVertices*(p-1);
            figure('name',pNames{p}); hold on;
            [im_srt,im_srt_idx]=sort(myInfoMapN(stOffset+1:stOffset+nVertices),'descend');
            isTopNVertex=false(nVertices,1);
            isTopNVertex(im_srt_idx(1:showTopNlabels))=true;
            for v=1:nVertices
                % get region coords from codebook
                cx=CB.center(v,1); cy=CB.center(v,2); cz=CB.center(v,3);
                thisImp=myInfoMapN(v+stOffset)*sphereScale;
                % scale spheres by their variable importance and translate them
                % to ROI center of mass
                sX=(thisImp*sXbase)+cx; sY=(thisImp*sYbase)+cy; sZ=(thisImp*sZbase)+cz;
                surf(sX,sY,sZ,'FaceAlpha',0.5,'CData',myCData,'edgecolor','none');
                if isTopNVertex(v)==true
                    text(cx,cy,cz,CB.sname{v}); % display abbreviation
                    %text(cx,cy,cz,CB.name{v}); % display full name
                end
            end
            % plot MNI space cortical mesh convex hull outline planes for
            % easier orientation
            for hi=1:2
                plot3(chop.Xs{hi}(chop.kas{hi}),chop.Ys{hi}(chop.kas{hi}),...
                    zeros(numel(chop.kas{hi}),1),'k--','LineWidth',2); % plot axial
                plot3(zeros(numel(chop.kss{hi}),1),chop.Ys{hi}(chop.kss{hi}),...
                    chop.Zs{hi}(chop.kss{hi}),'k--','LineWidth',2); % plot sagittal
                plot3(chop.Xs{hi}(chop.kcs{hi}),zeros(numel(chop.kcs{hi}),1),...
                    chop.Zs{hi}(chop.kcs{hi}),'k--','LineWidth',2); % plot coronal
            end
            % setup final props
            axis equal; axis tight; grid on; view(3); colormap gray;
            xlabel('left-right');
            ylabel('posterior-anterior');
            zlabel('ventral-dorsal');
            light;                      % make spheres shiny
        end
    otherwise
        error('Unknown approach');
end

%% ERROR ANALYSIS AND CLASSIFIER SELECTION
% after running several classifiers, we might want to check whether they
% yield the same results or not. Since the decisions are discrete and the
% test samples are the same, we need a paired testing procedure that
% accounts for this. We can use the McNemar test.
% 
% Important note: the fold structure should stay the same! This test only
% makes sense if the order of test samples is conserved between the test
% runs. So you may have to comment some of the code above.

% we can also compute the decision correctness vector for the classifier:
DR.(C_MYCLASSIFIER)=allPredictions==allTruths; % 1 means correct decision, 0 means error
% now we can run the CV code again and compare the new classifier to the
% old one, for instance
if isfield(DR,'NB') && isfield(DR,'SVM')
    [x2,isSignificantDiff]=jMcNemarTest(DR.NB,DR.SVM,0.05)
end