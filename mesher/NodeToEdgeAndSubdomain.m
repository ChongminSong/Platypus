function [ node2Edge, node2sd ] = NodeToEdgeAndSubdomain( meshEdge, edge2sd )
%Find edges and subdomains connected to each node

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

%ltx {\bf find subdomains connected to a node}
node2sd = cell(nNode,1); %ltx initilisation
for ii = 1:nNode %ltx loop over nodes
    node2sd{ii} = unique(vertcat(edge2sd{node2Edge{ii}}));
end

end

