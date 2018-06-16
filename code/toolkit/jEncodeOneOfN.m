function enc=jEncodeOneOfN(inp,varargin)
%
% encode a vector of data into a one-of-n encoding
%
% IN --     inp: a 1xT or Tx1 vector of integers
% OUT --    enc: a TxN matrix of one-of-N encoding
% 
% 2008 Nov - 0.1 - Jonas Richiardi
% - Initial release
% 2010 Apr JR
% - support for out-of-sample code values. This is useful if inp contains
% only a subset of the full code, but we want to have a full encoding
% anyways
% 0.2 2010 Oct JR
% - bug fix for non-contiguous alphabet case
% v0.3 2011 Nov JR
% - fix handling of case size(unique(inp)) ~= size(unique(fullAlphabet))

[fullAlphabet]=process_options(varargin,'fullAlphabet',[]);

inp=inp(:); % make column vector
fullAlphabet=fullAlphabet(:);
fullAlphaGen=[min(fullAlphabet):max(fullAlphabet)]';

%% remap if needed
if ~isempty(fullAlphabet)
    % check if the input has at least one representative of each member of the
    % alphabet set
    if all(size(unique(inp))==size(unique(fullAlphabet)))
        allAlphabetIsPresent=all(unique(inp)==unique(fullAlphabet));
    else
        allAlphabetIsPresent=false;
    end
    fullAlphabetIsContiguous=(length(fullAlphabet)==length(fullAlphaGen)) && (all(fullAlphabet==fullAlphaGen));
    if ~allAlphabetIsPresent
        % try to see if we need remapping        
        if (min(fullAlphabet)==1)
            if fullAlphabetIsContiguous
                % no need for remapping: alphabet is already contiguous and
                % starting with 1
                remapped=inp;
            else
                % need clever remapping
                error('Non-contiguous alphabet not supported when inp containing a subset of the full code');
            end
        else
            % need clever remapping
            error('Alphabet not starting with 1 is not supported when inp containing a subset of the full code');
        end
    else % all alphabet is present in input, do simple remapping
        [remapped,mapping]=remapDiscretesToContiguous(inp);
    end
else
    if length(inp)<=1
        warning('jEncodeOneOfN:smallInput',...
            'Input has length one or less and no fullAlphabet provided');
    end
    [remapped,mapping]=remapDiscretesToContiguous(inp);
end

%% do encoding
if ~isempty(fullAlphabet)
    %[remappedFA,mappingFA]=remapDiscretesForBNT(fullAlphabet(:));
    N=numel(unique(fullAlphabet));
else
    N=numel(unique(remapped));
end
T=numel(remapped);

enc=zeros(T,N);
for n=1:N
    enc(remapped==n,n)=1;
end
