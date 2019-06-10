function [coord, sdConn, sdSC] = combineTwoMeshes(d, onLine, res)
%cbgnltx 
%\myCodeTitle{Combine two non-matching meshes on a line defined by the function \texttt{onLine}}
%{\bf Inputs:}
%\begin{myDescription}
% \item[d\{i\}] mesh data of domain \texttt{i} 
% \item[onLine\{i\}]  handler of function definding the interface
% between domain \texttt{i} and \texttt{i+1} 
% \item[res] resolution of coordinates. Coordinates of the nodes on the 
% interface (except for the end nodes) will be rounded to mutiples
% of the resolution and nodes at the same coordinates will be merged.
% \end{myDescription}
%\vspace{1\baselineskip}\par
%{\bf Outputs:}
%\begin{myDescription}
%\item[coord(i,:)]  Coordinates of node \texttt{i}
%\item[sdConn\{isd,:\}(ie,:)] S-element conncetivity. The nodes of line
%   element \texttt{ie} in S-element \texttt{isd}.
%\item[sdSC(isd,:)] scaling centers of S-element \texttt{isd}.
% \end{myDescription}
%cendltx

%ltx find interface edges and nodes \label{Code_combineTwoMeshes_interface}
ctLine = cell(2,1); %ltx cells store the data of the 2 meshes at their interface
for ii = 1:2
    %ltx find interface edges (with both nodes on interface)
    ctLine{ii}.edgeID = intersect( ... 
          onLine(d{ii}.coord(d{ii}.meshEdge(:,1),:), res), ... 
          onLine(d{ii}.coord(d{ii}.meshEdge(:,2),:), res)); 
    edge = d{ii}.meshEdge(ctLine{ii}.edgeID,:)';
    %ltx round the coordinates of the nodes on interface to a specified tolerance
    %ltx to avoid inserting nodes very close to an existing nodes
    [node, ~, ic] = unique(edge); %ltx nodes on interface
    binCounter = histc(ic, 1:length(node)); 
    endNode = node(binCounter==1); %ltx nodes at the two ends of interface
    midNode = setdiff(node, endNode); %ltx middle nodes on interface
    d{ii}.coord(midNode,:) = round(d{ii}.coord(midNode,:)/res)*res;
    d{ii}.edgexy =  d{ii}.coord(edge,:);
    %ltx coordinates of nodes on interface
    ctLine{ii}.xy = d{ii}.coord(node,:); 
end

%cbgnltx
% \texttt{edgeXi\{1\}\{i\}} stores the parametric coordinates ($0< \xi < 1$) 
% of the points that will be inserted into the interface edge \texttt{i} of mesh 1 
% to match the interface nodes of mesh 2 (see Fig.~\ref{fig:ConnectingTwoRect}).
% Correspondingly, \texttt{edgeXi\{2\}\{i\}} is for
% the points to be inserted into the interface edges of mesh
% 2.\label{Code_combineTwoMeshes_insertNodes}
%cendltx
edgeXi{1} = findPointsOnLine(d{1}.edgexy, ctLine{2}.xy, res); 
edgeXi{2} = findPointsOnLine(d{2}.edgexy, ctLine{1}.xy, res);

%ltx subdivide the interface edges by inserting points to match the two meshes
for ii = 1:2
    edgeDiv = cell(length(d{ii}.meshEdge),1); %ltx initialisation
    %ltx parametric coordinates of points of subdivision at interface edges
    edgeDiv(ctLine{ii}.edgeID) = edgeXi{ii};
    n1 = length(d{ii}.coord)+1; %ltx number of nodes before subdivision
    %ltx subdivide edges and update the mesh of a domain
    [d{ii}.coord, d{ii}.sdConn ] = subdivideEdge(edgeDiv, ...
        d{ii}.coord, d{ii}.meshEdge, d{ii}.sdEdge);
    %ltx round the coordinates of newly created nodes on interface
    d{ii}.coord(n1+1:end,:) = round( ...
        d{ii}.coord(n1+1:end,:)/res )*res;
end

%ltx append the nodes of the two meshes \label{Code_combineTwoMeshes_update}
coord  = round([d{1}.coord; d{2}.coord]/(0.5*res)); 
[~, ia, ic] = unique(coord,'rows'); %ltx merge nodes
coord = coord(ia(:),:)*(0.5*res);
sdSC = [d{1}.sdSC; d{2}.sdSC];
%ltx append the S-element connectivity of meshes
sdConn = [ d{1}.sdConn; d{2}.sdConn];
%ltx shift nodal numbers in mesh 2 by the number of nodes in mesh 1 
shft = [ zeros(1,length(d{1}.sdConn)) ...
    length(d{1}.coord)*ones(1,length(d{2}.sdConn))]; 
%ltx update S-element connectivity
for ii = 1:length(sdConn)
    sdConn{ii} = ic(shft(ii)+sdConn{ii}); %ltx shifting
end

end