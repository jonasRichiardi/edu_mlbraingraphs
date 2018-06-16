function jAdjMatToNetworkFile(A,nNames,outFileName,varargin)
% dump an adjacency matrix to network file readable by cytoscape, Network
% workbench, Gephi, complex network package
%
% IN 
%   A: an nNodes x nNodes adjacency matrix for an undirected graph
%   nNames: node names, cell array of strings
%   outFileName: output file name
% 
% varargin
%     maxNodes: used for nNodes instead of size
%     formatOut: 'cytoScape', 'nwb', 'gephi', 'cnp', 'OSLOM'
% 
% Example use: jAdjMatToNetworkFile(A,num2strcell(1:5),'temp/testCynet.txt');
%   
% v1.0 Oct 2009 Jonas Richiardi
%   - support for cytoScape
% v1.1 Feb 2011 JR
%   - support for network workbench ( http://nwb.slis.indiana.edu )
% v1.2 Jul 2011 NL
%   - support for gephi ( http://gephi.org/ ), e.g. import into excel and
%   save as .csv file
%   (http://gephi.org/users/supported-graph-formats/spreadsheet/)
% v1.3 Mar 2013 JR
%   - support for Lev Muchnik's Complex Networks Package
%       http://www.levmuchnik.net/Content/Networks/ComplexNetworksPackage.html
% v1.3.1 May 2013 JR
% - support for OSLOM (http://ifisc.uib-csic.es/~jramasco/oslom/who.htm)
% v1.3.2 July 2013 JR
% - support for Cytoscape 3.0 (node attribute format changed)

% TODO: 
% vectorise all code as per OSLOM - huge speedup expected


% Note NWB can be used as interlingua (export to GraphML, Pajek .net/.mat,
% XGMML)

[maxNodes,formatOut]=process_options(varargin,'maxNodes',[],'formatOut','cytoScape');

%% sanity check
status=isSaneAdjacencyMatrix(A);
if ~isempty(status)
    error(status);
end

%% compute basic quantities
nNodes=size(A,1);
if ~isempty(maxNodes)
    nNodes=maxNodes;
end

degrees=sum(A,2);

switch formatOut
    case 'cytoScape'
        
        %% build edge definitions and strength

        outBufEdges{1}=['node1 node2 edgeVal'];
        for r=1:nNodes
            for c=r:nNodes
                if (A(r,c)~=0)
                    outBufEdges{end+1}=[nNames{r} ' ' nNames{c} ' ' num2str(A(r,c)) ]; % for Gephi add:  ' ' ' "Undirected"'
                end
            end
        end
        
        %% build node attributes file
        outBufNodes=cell(nNodes+1,1);
        outBufNodes{1}='nodeVal';
        for n=1:nNodes
            outBufNodes{n+1}=[nNames{n} ' = ' num2str(degrees(n))];
        end
        
        %% dump edge attributes file
        fid=fopen(outFileName,'w');
        for i=1:numel(outBufEdges)
            fprintf(fid,'%s\n',outBufEdges{i});
        end
        fclose(fid);
        
        %% dump node attributes file (Cytoscape v 2.x)
        fid=fopen([outFileName '.noa'],'w');
        for i=1:numel(outBufNodes)
            fprintf(fid,'%s\n',outBufNodes{i});
        end
        fclose(fid);
        
        %% dump node attributes file (Cytoscape v 3.0.x)
        fid=fopen([outFileName '.noa3.txt'],'w');
        outBufNodes{1}='nodeName nodeVal';
        for i=1:numel(outBufNodes)
            % replace equal char on the fly with space char
            outBufNodes{i}=strrep(outBufNodes{i},'= ','');
            fprintf(fid,'%s\n',outBufNodes{i});
        end
        fclose(fid);
        
    case 'nwb'
        outBuf=cell(0);
        % buffer nodes
        outBuf{1}=['*Nodes ' num2str(nNodes)];
        outBuf{2}=['id*int label*string'];
        for n=1:nNodes
            outBuf{end+1}=[num2str(n) ' "' nNames{n} '"'];
        end
        % buffer edges
        outBuf{end+1}=['*UndirectedEdges ' num2str(nchoosek(nNodes,2))];
        outBuf{end+1}=['source*int target*int corrWeight*float'];
        % inefficient! could do in three steps of vectoring
        for l=1:nNodes
            for c=(l+1):nNodes
                theWeight=round(A(l,c)*100);
                outBuf{end+1}=[num2str(l) ' ' num2str(c) ' ' num2str(theWeight)];
            end
        end
        % dump to file
        fid=fopen([outFileName '.nwb'],'w');
        for i=1:numel(outBuf)
            fprintf(fid,'%s\n',outBuf{i});
        end
        fclose(fid);
        
    case 'gephi'
        
        % build edge definitions and strength

        outBufEdges{1}=['Source Target Weight Type'];
        for r=1:nNodes
            for c=r:nNodes
                if (A(r,c)~=0)
                    outBufEdges{end+1}=[nNames{r} ' ' nNames{c} ' ' num2str(A(r,c)) ' ' ' "Undirected"'];  
                end
            end
        end
        
        % build node attributes file
%         outBufNodes=cell(nNodes+1,1);
%         outBufNodes{1}='nodeDegrees';
%         for n=1:nNodes
%             outBufNodes{n+1}=[nNames{n} ' = ' num2str(degrees(n))];
%         end
%         
        % dump edge attributes file
        fid=fopen(outFileName,'w');
        for i=1:numel(outBufEdges)
            fprintf(fid,'%s\n',outBufEdges{i});
        end
        fclose(fid);
        
        % dump node attributes file
%         fid=fopen([outFileName '.noa'],'w');
%         for i=1:numel(outBufNodes)
%             fprintf(fid,'%s\n',outBufNodes{i});
%         end
%         fclose(fid);
    case 'cnp'
        % wants 1 row and col as headers
        outBuf=[[0 1:nNodes]' [1:nNodes ; A]];
        save(outFileName,'outBuf','-ascii');
        
    case 'OSLOM'
        zeroDeg_idx=sum(A)==0;  % find zero-degree vertices
        NZ=A~=0;                % compute logical mask of non-zero edge labels
        NZidx=find(NZ);         % ... and its corresponding index
        NZsub=ind2subv(size(A),NZidx); % ... and the corresponding subscripts
        % get non-empty pairs
        names_r=nNames(NZsub(:,1))'; % name for vertices in rows of adjmat
        names_c=nNames(NZsub(:,2))'; % same for columns
        vals=num2strcell(full(A(NZidx)))';  % edge label
        % add zero-degree vertices. They need to be output with
        % a tiny self-weight so that at least their IDs is part of the vertex
        % set and OSLOM does not complain that the vertex IDs in the hint
        % file (class labels) do not match those in the adjmat file.
        % see OSLOM's undirected_network.h : static_network::translate
        names_r=[names_r;nNames(zeroDeg_idx)'];
        names_c=[names_c;nNames(zeroDeg_idx)'];
        vals=[vals;num2strcell(repmat(eps,sum(zeroDeg_idx),1))'];
        spaces=repmat({' '},numel(names_r),1);
        outBuf=strcat(names_r,spaces,names_c,spaces,vals);
        
        %% dump edge file
        fid=fopen(outFileName,'w');
        for i=1:numel(outBuf)
            fprintf(fid,'%s\n',outBuf{i});
        end
        fclose(fid);
        
    otherwise
        error('Unknown output format');
end