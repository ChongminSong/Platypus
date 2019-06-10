function [ meshEdge, sdEdge] = findMeshEdge( sdConn )
%Construct mesh connectivity data
%
%Input
%   sdConn: subdomain/element connectivity
%      when sdConn is a matrix, 
%           sdConn(i,:) are the nodes of element i  
%      when sdConn is a cell array of a 1D array, 
%           sdConn{i} are the nodes of ploygon i
%      when sdConn is a cell array of a matrix, 
%           sdConn{i} are the nodes of line elements of a subdomain
%           
%Output
%   meshEdge(i,1:2)  - the 2 nodes on edge i in a mesh.
%   sdEdge{i}        - edges forming subdomain i
%   edge2sd{i}       - subdomains connected to edge i
%   node2Edge{i}     - edges connected to node i
%   node2sd{i}       - subdomains connected to node i

nsd = length(sdConn); %number of subdomains/elements

%ltx {\bf construct connectivity of edges to subdomains/elements}
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
        i1 = i2 + 1; %ltx update the counter for the next element
    end
    %ltx combine edges of all elements. The 2 nodes of an edge are sorted
    %ltx in ascending order. Each edge is stored as one row.
    meshEdge = (sort([meshEdge{:}]))';
else    %ltx input of \texttt{sdConn} is a cell array
    sdEdge = cell(nsd,1); %ltx initialisation
    i1 = 1;
    if size(sdConn{1},1) > 1 
        %ltx SBFEM subdomains (a cell has multiple rows)
        %ltx all the edges of subdomains are numbered in this loop
        for ii = 1:nsd %ltx loop over subdomains
            i2 = i1 + length(sdConn{ii}) - 1; %ltx count the edges 
            sdEdge{ii} = i1:i2; %ltx edges of a subdomain
            i1 = i2 + 1;
        end
         %ltx combine edges of all polygons. The 2 nodes of an edge are 
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
            i1 = i2 + 1; %ltx update the counter for the next element
        end
        %ltx combine all edges of all elements. The 2 nodes of an edge are
        %ltx sorted in ascending order. Each edge is stored as one row.
        meshEdge = (sort([meshEdge{:}]))';
    end
end

%ltx remove duplicated entries of edges
[meshEdge, ~, ic] = unique(meshEdge,'rows');
for ii = 1:nsd %ltx loop over subdomains/elements
    sdEdge{ii} = ic(sdEdge{ii}); %ltx update edge numbers
end

end
