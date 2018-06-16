function [x2,isSignificant]=jMcNemarTest(DR1,DR2,alpha)
% [x2,isSignificant]=jMcNemarTest(DR1,DR2,alpha)
% 
% performs the McNemar hypothesis test on the significance of differences
% between two 2-class classifiers, as per Kuncheva "combining
% classifiers" book p.14.
%
% IN - DR1: nTests x 1 column vector of {0,1} flags corresponding to
%           {wrong, correct} classifier 1 decisions
%      DR2: likewise for classifier 2
%      alpha: 1-significance level of the test (typ 0.05)
% OUT - x2: represents the discrepency between the expected counts of
%           errors if the classifiers perform the same (i.e., both have
%           the same number of errors, which is the null hypothesis), and
%           the actual counts
%       isSignificant: true if the difference between the 2 classifiers is
%           significant at the significance level, false otherwise
%
% v1.0 Nov 2006 Jonas Richiardi
% v1.1 Dec 2006 JR for MCS/thesis
% - check for zero diff before division
% v1.1.1 May 2010 JR 
% - code and doc cleanup for advanced fMRI course

if ~((alpha<=1) && (alpha>=0))
    error('jMcNemarTest: alpha must be between 0 and 1')
end

% make column vectors
DR1=DR1(:);
DR2=DR2(:);

% compute critical value
critical=chi2inv(1-alpha,1);

% compute number of times classifier 1 is wrong while classifier 2 is
% correct
N01=sum(~DR1&DR2);
% compute number of times classifier 2 is wrong while classifier 1 is
% correct
N10=sum(DR1&~DR2);

% compute discrepancy between expected counts for null hypothesis
% and observed counts
if ((N01==0) && (N10==0))
    x2=0;
else
    x2=((abs(N01-N10)-1)^2)/(N01+N10);
end

if (x2 > critical)
    isSignificant=true;
else
    isSignificant=false;
end
