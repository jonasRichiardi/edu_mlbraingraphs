function bf=graphSelectFeatures(X,cLabs,D,checkWithinClassVar,doPlot)
% function bf=graphSelectFeatures(X,cLabs,D,checkWithinClassVar,doPlot)
% 
% Select best features using point-biserial correlation criterion and
% univariate ranking search method.
%
% This is not necessarily a good idea in practice, but serves as an
% illustration.
%
% IN
%   X:      nFeats x nPatterns matrix of column feature vectors
%   cLabs:  nPatterns x 1 vector of class labels
%   D:      scalar indicating how many features to select
%   checkWithinClassVar: boolean, whether to check if within-class variance
%           is not pathological
%   doPlot: boolean, whether to plot sorted features graph or not
%
% OUT
%   bf:     D x 1 vector of feature indices indicating the best features
%
% v1.0 July 2013 Jonas Richiardi for Grenoble summer school

% sanity check
[nFeats, nPats]=size(X);
if length(cLabs)~=nPats
    error('cLabs length inconsistent with number of patterns');
end
if D>nFeats
    error('D must be smaller or equal to number of features in X');
end

c1idx=cLabs==0;
c2idx=cLabs==1;

% compute per-feature point-biserial correlation (which is a special
% case of pearson product-moment correlation, so no need to implement
% your own code except for maybe speed gains).
% Use for-loop to save memory over the matrix-matrix pairwise corr case
Rpb=zeros(1,nFeats);
parfor f=1:nFeats
    % check pathological cases of all features with same value across
    % dataset
    if (std(X(f,:))==0)
        Rpb(f)=-inf;
    % if needed, check if std is zero across a particular class (probably indicates
    % faulty feature extraction)
    elseif checkWithinClassVar==true && ((std(X(f,c1idx))==0 || std(X(f,c2idx))==0))
        Rpb(f)=-inf;
    % not a pathological case, do computation
    else
        rpb=corrcoef(X(f,:),cLabs);
        Rpb(f)=abs(rpb(1,2)); % extract correlation coeff from matrix
    end
end

% get ranking
[Rpb_srt,Rpb_srt_idx]=sort(Rpb,'descend');

% extract top D feature indices
bf=Rpb_srt_idx(1:D);

% check if we have some pathological data still
isPathological=isinf(Rpb_srt(1:D));
if any(isPathological)
    % warn the user to deal with it by reducing
    % we could also do it here automatically but maybe not desirable
    warning('Some features are still pathological, try reducing dimension further');
    disp(['Look at features ' mat2str(find(isPathological))]);
end

if doPlot==true
    figure; plot(Rpb_srt,'LineWidth',2.5);
    xlabel('Feature number');
    ylabel('|r_{pb}|'); ylim([0 1]);
    grid on; hold on;
    % D-line
    line([D D],[0 1],'LineStyle','--','LineWidth',2,'Color','r');
end