function [ meshEdge, sdEdge, edge2sd, node2Edge, node2sd] = ...
    meshConnectivity( sdConn )
%Construct mesh connectivity data
%
%Input
%   sdConn: S-element/element connectivity
%      when sdConn is a matrix, 
%           sdConn(i,:) are the nodes of element i  
%      when sdConn is a cell array of a 1D array, 
%           sdConn{i} are the nodes of ploygon i
%      when sdConn is a cell array of a matrix, 
%           sdConn{i} are the nodes of line elements of S-element i
%           
%Output
%   meshEdge(i,1:2)  - the 2 nodes on line i in a mesh.  
%                      The node number of the first node is 
%                      larger than that of the 2nd one
%   sdEdge{i}        - lines forming S-element i, 
%                      >0 when an line follows anti-clockwise 
%                         direction around the scaling centre.   
%                      <0 otherwise 
%   edge2sd{i}       - S-elements connected to edge i
%   node2Edge{i}     - edges connected to node i
%   node2sd{i}       - S-elements connected to node i

nsd = length(sdConn); %ltx number of S-elements/elements

%ltx {\bf construct connectivity of edges to S-elements/elements}
%ltx \texttt{sdConn} is intepreted based on its data type
if ~iscell(sdConn) %ltx \texttt{sdConn} is a matrix (triangular elements) 
    %ltx the following loop collects the edges of all elements
    sdEdge = cell(nsd,1); %ltx initialisation
    i1 = 1; %ltx counter of edges
    meshEdge = cell(nsd,1); %ltx initialisation
    for ii = 1:nsd %ltx loop over elements
        eNode = sdConn(ii,:); %ltx nodes of an element
        i2 = i1 + length(eNode) - 1; %ltx count the edges 
        %ltx element edges. Each edge is defined by 2 nodes on a column
        meshEdge{ii} = [eNode; eNode(2:end) eNode(1)]; 
        sdEdge{ii} = i1:i2; %ltx edges of an element
        %ltx store the edges to be reversed in sorting as negative values
        idx = find( meshEdge{ii}(1,:) > meshEdge{ii}(2,:) );
        sdEdge{ii}(idx) = - sdEdge{ii}(idx); 
        i1 = i2 + 1; %ltx update the counter for the next element
    end
    %ltx combine edges of all elements. The 2 nodes of an edge are sorted
    %ltx in ascending order. Each edge is stored as one row.
    meshEdge = (sort([meshEdge{:}]))';
else    %ltx input of \texttt{sdConn} is a cell array
    sdEdge = cell(nsd,1); %ltx initialisation
    i1 = 1;
    if size(sdConn{1},1) > 1 
        %ltx S-elements (a cell has multiple rows)
        %ltx all the edges of S-elements are numbered in this loop
        for ii = 1:nsd %ltx loop over S-elements
            i2 = i1 + length(sdConn{ii}) - 1; %ltx count the edges 
            sdEdge{ii} = i1:i2; %ltx edges of an S-element
            %ltx store the edges to be reversed in sorting as negative values
            %ltx each element edge is defined by 2 nodes as a row
            idx = find( sdConn{ii}(:,1) > sdConn{ii}(:,2) );
            sdEdge{ii}(idx) = - sdEdge{ii}(idx); 
            i1 = i2 + 1;
        end
         %ltx combine edges of all S-elements. The 2 nodes of an edge are 
         %ltx sorted in ascending order. Each edge is stored as one row.
        meshEdge = sort(vertcat(sdConn{:}),2);
    else
        %ltx polygon element (closed loop specified by vertices)
        %ltx the following loop collects the edges of all polygons
        meshEdge = cell(nsd,1); %ltx initialisation
        for ii = 1:nsd %ltx loop over polygons
            eNode = sdConn{ii}; %ltx nodes of a polygon
            i2 = i1 + length(eNode) - 1; %ltx count the edges 
            %ltx each element edge is defined by 2 nodes as a column
            meshEdge{ii} = [eNode; eNode(2:end) eNode(1)]; 
            sdEdge{ii} = i1:i2; %ltx edges of a polygon
            idx = find( meshEdge{ii}(1,:) > meshEdge{ii}(2,:) );
            %ltx edge to be reversed
            sdEdge{ii}(idx) = -sdEdge{ii}(idx); 
            i1 = i2 + 1; %ltx update the counter for the next element
        end
        %ltx combine all edges of all elements. The 2 nodes of an edge are
        %ltx sorted in ascending order. Each edge is stored as one row.
        meshEdge = (sort([meshEdge{:}]))';
    end
end

%ltx remove duplicated entries of edges
[meshEdge, ~, ic] = unique(meshEdge,'rows');
for ii = 1:nsd %ltx loop over S-elements/elements
    %ltx update edge numbers
    sdEdge{ii} = sign(sdEdge{ii}(:)).*ic(abs(sdEdge{ii})); 
end

if nargout < 3 
    return;
end

%ltx {\bf find S-elements/elements connected to an edge}
a = abs(cell2mat(sdEdge)); %ltx edges of all S-elements/elements
%ltx the following loop matchs S-element/element numbers to edges
asd = zeros(length(a),1); %ltx initialisation
ib = 1; %ltx pointer to the first edge of an S-element/element
for ii = 1:nsd %ltx loop over S-elements/elements
    ie = ib + length(sdEdge{ii}) - 1;
    asd(ib:ie) = ii; %ltx  edge a(i) is connected to S-element/element asd(i) 
    ib = ie + 1; %update the pointer
end
%ltx sort S-element numbers according to edge number
[c, indx] = sort(a); asd = asd(indx); 
%ltx the following loop collects the S-elements connected to nodes
ib = 1; %ltx pointer to the 1st S-element/element connected to an edge
nMeshedges = length(meshEdge); %ltx number of edges in mesh
edge2sd = cell(nMeshedges,1); %ltx initilisation 
for ii = 1:nMeshedges-1 %ltx loop over edges  (except for the last one)
    if c(ib+1) == ii %ltx two S-elements/elements connected to an edge
        edge2sd{ii} = asd(ib:ib+1); %ltx store the S-elements/elements 
        ib = ib + 2; %ltx update the pointer
    else %ltx one S-element/element connected to an edge
        edge2sd{ii} = asd(ib); %ltx store the S-element/element
        ib = ib + 1; %ltx update the pointer
    end
end
%ltx the S-elements/elements connected to the last edges
edge2sd{ii+1} = asd(ib:end);

if nargout < 4 
    return;
end

%ltx {\bf find edges connected to a node}
a = reshape(meshEdge',1,[]); %ltx nodes on edges
%ltx edge numbers correponding to nodes in \texttt{a}
edgei = reshape([1:size(meshEdge,1); 1:size(meshEdge,1)],1,[]); 
%ltx sort edge number according to node number
[c, indx] = sort(a); edgei = edgei(indx);
ib = 1; %ltx pointer to the 1st edge connected to a node
nNode = c(end); %ltx number of nodes
node2Edge = cell(nNode,1); %ltx initilisation
for ii = 1:nNode-1 %ltx loop over nodes (except for the last one)
    %ltx pointer to the last edge connected to a node
    ie = ib-2+find(c(ib:end)~=ii, 1, 'first'); 
    node2Edge{ii} = edgei(ib:ie); %ltx store edges connected to node
    ib = ie + 1; %ltx update the pointer
end
%ltx store edges connected to the last node
node2Edge{nNode} = edgei(ib:end); 

%ltx {\bf find S-elements connected to a node}
node2sd = cell(nNode,1); %ltx initilisation
for ii = 1:nNode %ltx loop over nodes
    node2sd{ii} = unique(vertcat(edge2sd{node2Edge{ii}}));
end

end
