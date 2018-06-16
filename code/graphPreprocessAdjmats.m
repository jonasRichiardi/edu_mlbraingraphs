function CM3Dout=graphPreprocessAdjmats(CM3D,mod,classNames,threshMethod,threshVal)
% function CM3Dout=graphPreprocessAdjmats(CM23D,mod,threshMethod,threshVal)
% preprocess adjacency matrices by removing self-loops, checking they
% are otherwise legal, and thresholding them
%
% IN
%
% OUT
%
% VERSION
% v1.0 July 2013 Jonas Richiardi for Grenoble Summer School
% - initial release

%forceSymmetry=true; % shoudl we force symmetry if original adjmat is non-symmetric?

% remove self-loops, and then check in more detail for adjacency matrix
% validity, fix asymmetry if present
CM3Dout=struct;
CM3Dout.(mod)=struct;
for c=classNames
    nAdjmats=size(CM3D.(mod).(c{:}),3);
    if size(CM3D.(mod).(c{:}),1)~=size(CM3D.(mod).(c{:}),2)
        error('adjmat must be square');
    end
    nVertices=size(CM3D.(mod).(c{:}),1);
    CM3Dout.(mod).(c{:})=zeros(nVertices,nVertices,nAdjmats);
    for a=1:nAdjmats
        % remove diagonal
        tmpMat=CM3D.(mod).(c{:})(:,:,a)-diag(diag(CM3D.(mod).(c{:})(:,:,a)));
        % take care of possible nans
        tmpMat(isnan(tmpMat))=0;
        % ready to check sanity now
        checkResult=isSaneAdjacencyMatrix(tmpMat);
        if ~isempty(checkResult)
            disp(checkResult);
            if ~isempty(strfind(checkResult,'SYMMETRIC'))
                disp(['Attempting to fix adjmat ' c{:} '/' num2str(a) ...
                    ' by symmetrising...']);
                tmpMat=(tmpMat+tmpMat.')/2;
                checkResult2=isSaneAdjacencyMatrix(tmpMat);
                if ~isempty(checkResult2)
                    error(['Adjmat ' c{:} '/' num2str(a) ...
                        ' is still not valid. Giving up.']);
                else
                    disp('... fixing worked.');
                    % seems OK now, copy to output
                    CM3Dout.(mod).(c{:})(:,:,a)=tmpMat;
                end
            else
                error(['Adjmat ' c{:} '/' num2str(a) ...
                    ' is not valid. Giving up.']);
            end
        else
            % seems OK, copy to output
            CM3Dout.(mod).(c{:})(:,:,a)=tmpMat;
        end
    end
end

% now we can threshold
switch threshMethod
    case ''
        disp('No thresholding applied');
    case 'edgeDensity'
        for c=classNames
            nAdjmats=size(CM3D.(mod).(c{:}),3);
            nVertices=size(CM3D.(mod).(c{:}),1);
            for a=1:nAdjmats
                thisA=CM3Dout.(mod).(c{:})(:,:,a);
                allVals=jUpperTriMatToVec(thisA,1);
                myPct=prctile(allVals,[100-threshVal*100 threshVal*100]);
                CM3Dout.(mod).(c{:})(:,:,a)=thisA.*double(thisA>=myPct(1));
                EDachieved=sum(sum(CM3Dout.(mod).(c{:})(:,:,a)~=0))/numel(thisA);
                fprintf('Edge density target: %1.2f, achieved: %1.2f\n',...
                    threshVal,EDachieved);
                if abs(threshVal-EDachieved)>=0.2
                    warning('GPA:EDtargetMissed',['Edge density target '...
                        ' could not be met for ' c{:} '/' num2str(a) ]);
                    fprintf('The unthresholded edge density is %1.2f, try setting a lower target.\n',nnz(thisA)/numel(thisA));
                    fprintf('%2.1f percent of edge labels are unique\n',100*numel(unique(allVals))/numel(allVals));
                end
            end
        end
    otherwise
        error(['Unknown thresholding method ' threshMethod]);
end
