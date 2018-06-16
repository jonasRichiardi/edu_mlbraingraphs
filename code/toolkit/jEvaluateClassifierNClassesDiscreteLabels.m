function resultStats=jEvaluateClassifierNClassesDiscreteLabels(truth,hatv,varargin)
% function resultStats=jEvaluateClassiferNCLasses(truth,hatv,varargin)
% give accuracy figures etc on the classifier's output for N-class problems
% uses netlab
% IN:   truth - a nCases x 1 vector of correct class labels
%       hatv - a nCases x 1 vector of classified class labels
% OUT:  stats about the classifier's performance
%       . correctRate: correct recognition rate
%       . confMat: confusion matrix
% v1.0 oct 2009 Jonas Richiardi
% - initial release
% v1.1 Oct 2010 JR
% - bug fix for non-contiguous classes case

truth=truth(:);
hatv=hatv(:);

nCases=size(truth,1);

[plotStuff,nClasses,fullAlphabet]=process_options(varargin,'plotStuff',0,...
    'nClasses',[],'fullAlphabet',[]);

% guess fullAlphabet from nClasses
if isempty(fullAlphabet) && isempty(nClasses)
    if nCases==1
        error('Cannot guess how many classes there are');
    else
        nClasses=numel(unique(truth)); % assume all classes are represented in the test set
        warning('jecncd:guessing',...
            ['Guessing there are ' num2str(nClasses) ' classes.']);
    end
elseif isempty(fullAlphabet) && ~isempty(nClasses)
    fullAlphabet=1:nClasses;
elseif ~isempty(fullAlphabet) && isempty(nClasses)
    nClasses=numel(unique(fullAlphabet));
    if nClasses~=numel(fullAlphabet)
        error(['Bad fullAlphabet specification: there are repeating ' ...
            'elements']);
    end
elseif ~isempty(fullAlphabet) && ~isempty(nClasses)
    if nClasses~=numel(unique(fullAlphabet))
        error(['nClasses is not consistent with fullAlphabet']);
    end
end
    

truth1OfN=jEncodeOneOfN(truth,'FullAlphabet',fullAlphabet);
hatv1OfN=jEncodeOneOfN(hatv,'FullAlphabet',fullAlphabet);

% use netlab to compute confmat
[myConfmat,correctRate]=confmat(hatv1OfN,truth1OfN); % true x pred
if (plotStuff==1)
    conffig(hatv1OfN,truth1OfN);
end
% compute per-class accuracy and total errors
cAcc=zeros(1,nClasses);
for c=1:nClasses
    cAcc(c)=myConfmat(c,c)/sum(myConfmat(c,:));
end

resultStats.correctRate=correctRate(1);
resultStats.confMat=myConfmat;
resultStats.classAccs=cAcc;
resultStats.meanClassAcc=mean(cAcc);
resultStats.nErrors=sum(sum(myConfmat))-sum(diag(myConfmat)); % errors are off-diagonal terms

% compute positive and negative predctive values if only two classes
% assuming by convention that first class is healthy and second is patient
if nClasses==2
    TN=myConfmat(1,1); FP=myConfmat(1,2);
    FN=myConfmat(2,1); TP=myConfmat(2,2);
    resultStats.PPV=TP/(TP+FP);     % sens: TP/(TP+FN)
    resultStats.NPV=TN/(TN+FN);     % spec: TN/(TN+FP)
end
