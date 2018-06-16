function [out,mapping]=remapDiscretesToContiguous(myData)
% remap disjoint discrete values into contiguous discrete values
% starting from 1.
%
% IN -- myData: a TxN array of row vectors. Each column is a random
%           variable.
% OUT -- out: a TxN array of remapped row vectors.
% v1.0 May 2008 Jonas Richiardi

if(size(myData,2)>1)
    error('Code not ready for multivariate data yet');
end

baseAlphabet=unique(myData);
remapAlphabet=(1:length(baseAlphabet))';

out=myData;
% ugly loop solution
for s=1:length(baseAlphabet)
    out(myData==baseAlphabet(s))=remapAlphabet(s);
end

mapping=[baseAlphabet remapAlphabet];