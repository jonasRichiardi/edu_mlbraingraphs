function imps=jNBFeatureImportance(TR,cLabels,varargin)
% compute importance of feature for naive Bayes classifier
% first approximation: importance is related to log-odds ratio a feature
% (if using normal densities), or its distribution-free Equal Error Rate
% (if using a less-parametric density estimation method)
% 
% IN 
%   TR : nFeats x nSamples matrix of column vectors
%   cLabels: nSamples x 1 vector of class labels (1,2)
% v0.1 beta Dec 2011 Jonas Richiardi

[method]=process_options(varargin,'method','oddsRatio');

cLabels=cLabels(:);
if ~all(unique(cLabels)==[1 2]')
    error('labels must be 1 and 2. Only supports two-class probs.')
end

nFeats=size(TR,1);
nSamples=size(TR,2);
imps=zeros(1,nFeats,'single');
parfor d=1:nFeats
    if mod(d,500)==0
        fprintf('%s ', num2str(d));
    end
    switch method
        case 'EER'
            error('EER disabled in this version');
            [EER, EERthreshold]=jEER_DET(TR(d,cLabels==1),TR(d,cLabels==2),...
                'plotDETcurve',0,'nSteps',101,'verbose',0);
            imps(d)=EER;
        case 'oddsRatio'
            [c1mu, c1sigma] = normfit(TR(d,cLabels==1));
            [c2mu, c2sigma] = normfit(TR(d,cLabels==2));
            lor=0;
            for s=1:nSamples
                c1LL=-normlike([c1mu c1sigma],TR(d,s));
                c2LL=-normlike([c2mu c2sigma],TR(d,s));
                lor=lor+abs(c1LL-c2LL);
            end
            lor=lor/nSamples;   % normalise
            imps(d)=lor;
        otherwise
            error('unknown method')
    end
end