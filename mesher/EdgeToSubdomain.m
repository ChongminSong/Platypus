function [ edge2sd ] = findEdgeToSubdomainConnetivity( sdEdge )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%ltx {\bf find subdomains/elements connected to an edge}
a = cell2mat(sdEdge); %ltx accumulate edges of all subdomains/elements
%ltx the following loop matchs subdomain/element numbers to edges
asd = zeros(length(a),1); %ltx initialisation
ib = 1; %ltx pointer to the first edge of a subdomain/element
for ii = 1:nsd %ltx loop over subdomains/elements
    ie = ib + length(sdEdge{ii}) - 1;
    asd(ib:ie) = ii; %ltx  edge a(i) is connected to subdomain/element asd(i) 
    ib = ie + 1; %update the pointer
end
%ltx sort subdomain numbers according to edge number
[c, indx] = sort(a); asd = asd(indx); 
%ltx the following loop collects the subdomains connected to nodes
ib = 1; %ltx pointer to the 1st subdomain/element connected to an edge
nMeshedges = length(meshEdge); %ltx number of edges in mesh
edge2sd = cell(nMeshedges,1); %ltx initilisation 
for ii = 1:nMeshedges-1 %ltx loop over edges  (except for the last one)
    if c(ib+1) == ii %ltx 2 subdomains/elements connected to an edge
        edge2sd{ii} = asd(ib:ib+1); %ltx store the subdomains/elements 
        ib = ib + 2; %ltx update the pointer
    else %ltx 1 subdomain/element connected to an edge
        edge2sd{ii} = asd(ib); %ltx store the subdomain/element
        ib = ib + 1; %ltx update the pointer
    end
end
%ltx the subdomains/elements connected to the last edges
edge2sd{ii+1} = asd(ib:end);

end

