function myMap=jPlotDecisionBoundary(myC,cType,TR)
% function myMap=jPlotDecisionBoundary(myC,cType,TR)
% plot decision boundary for a classifier with 2 features
%
% IN
%   myC: a classifier built by decoding_practical.m
%   cTYpe: string, classifier type ('NB','SVM','RF');
%   TR: training set (column vectors of features)
% OUT
%   myMap.MLDecisions: a map of classifier decisions per grid point
%
% v1.0 2008 Jonas Richiardi
% - initial release
% v1.1 May 2010 Jonas Richiardi for advanced fMRI course
% - cleanup for course

%% synthesize data
minVal1=min(TR(1,:)); % min val for feature 1
maxVal1=max(TR(1,:)); % max val for feature 1
minVal2=min(TR(2,:)); % min val for feature 2
maxVal2=max(TR(2,:)); % max val for feature 2

minVal=min(minVal1,minVal2);
maxVal=max(maxVal1,maxVal2);

% spread data on a grid
plotMargin=0.2;
nGridPoints=40;
x = linspace(minVal1-abs(minVal1*plotMargin),maxVal1+abs(maxVal1*plotMargin),nGridPoints);
y = linspace(minVal2-abs(minVal2*plotMargin),maxVal2+abs(maxVal2*plotMargin),nGridPoints);
[myMap.XX, myMap.YY] = meshgrid(x,y);


%% compute decision at each grid point
% ugly slow loop!
disp('Computing decisions over input feature space...');
myMap.MLDecisions=zeros(size(myMap.XX));
for ri=1:size(myMap.XX,1)
    for ci=1:size(myMap.YY,2)
        tmpFeatVec=[myMap.XX(ri,ci);myMap.YY(ri,ci)];
        switch cType
            case 'NB'
                myMap.MLDecisions(ri,ci)=predict(myC,tmpFeatVec');
            case 'SVM'
                myMap.MLDecisions(ri,ci)=svmpredict(0,tmpFeatVec',myC);
            case 'DT'
                [predLabels,~] = myC.eval(tmpFeatVec');
                myMap.MLDecisions(ri,ci)=cell2num(predLabels);
            case 'RF'
                myMap.MLDecisions(ri,ci) = classRF_predict(tmpFeatVec',myC);
            otherwise
                error('Unknown classifier');
        end
    end
end

%% plot decision boundary
hold on;
% Compute contour of decision regions
[C,h]=contour(myMap.XX,myMap.YY,double(myMap.MLDecisions),1);
set(h,'LineWidth',3.0,'LineStyle',':','Color',[0 0 0]);
