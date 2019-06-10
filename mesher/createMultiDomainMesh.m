function [coord,sdConn,sdSC] = createMultiDomainMesh(dDefs, para)
%cbgnltx 
%\myCodeTitle{% 
%Create a mesh of multiple domains by combining meshes 
% generated independently for each domain}
%{\bf Inputs:}
%\begin{myDescription}
%\item[dDefs] structure variable storing the data of the problem domains
%and the interfaces in the following fields: \\
%\begin{myTabular}{lll}% 
%\setlength{\itemsep}{0pt}%
%  \texttt{prodef\{i\}}: & handle of problem definition function of domain \texttt{i}\\
%  \texttt{cshft(i,:)}: & amount of translation of coordinates of domain \texttt{i}\\
%  \texttt{onLine\{i\}}: & handle of function defining the interface between\\
%                         &  domain \texttt{i-1} and domain \texttt{i}  \\
%\end{myTabular}
%\item[para\{i\}] cell array of structure variables containing parameters 
%      for mesh generation of domain \texttt{i}. The fields are:\\
%\begin{myTabular}{lll}
%   \texttt{mesher}: &   = 1,& DistMesh \\
%          &  = 2,& PolyMesher \\
%          &  = 3,& direct input of triangular mesh\\
%          &  = 4,& direct input of polygon mesh \\
%          &  = 'keyword',& access user-defined function to input \\
%          &       &or generate S-element mesh\\
%   \texttt{h0}:& & initial edge length (DistMesh only) \\
%   \texttt{nPolygon}:& & number of polygons (PolyMesher only)\\
%   \texttt{nDiv}:& & number of 2-node line elements per edge
% \end{myTabular}
%\end{myDescription}
%\vspace{1\baselineskip}\par
%{\bf Outputs:}
%\begin{myDescription}
%\item[coord(i,:)]  Coordinates of node \texttt{i}
%\item[sdConn\{isd,:\}(ie,:)] S-element conncetivity. The nodes of line
%   element \texttt{ie} in S-element \texttt{isd}.
%\item[sdSC(isd,:)] scaling centers of S-element \texttt{isd}.
%\end{myDescription}
%cendltx

nDomain = length(dDefs.prodef); %ltx number of domains 
d = cell(nDomain,1); %ltx initilise cell variable storing domains
minLeng = 1d20; %ltx initilise minimum length of edges in a mesh
for id = 1:nDomain %ltx loop over domain \label{createMultiDomainMeshDomain}
    [ coord, sdConn, sdSC ] = createSBFEMesh( ...
        dDefs.prodef{id}, para{id}.mesher, para{id});
    [ meshEdge, sdEdge] = meshConnectivity( sdConn );
    %ltx shift coordinates of nodes and scaling centres
    coord = bsxfun(@plus, coord, dDefs.cshft(id,:));
    sdSC = bsxfun(@plus, sdSC, dDefs.cshft(id,:));
    d{id} = struct('coord',coord, 'sdConn',{sdConn}, 'sdSC',sdSC, ...
         'meshEdge',{meshEdge}, 'sdEdge',{sdEdge});
    %ltx length of shortest edge in mesh
    dxy = coord(meshEdge(:,1),:) - coord(meshEdge(:,2),:);
    minLeng = min( minLeng , sqrt(min(sum(dxy.*dxy,2))) );
% %     figure; % plot mesh
% %     opt=struct('sdSC', d{id}.sdSC, 'PlotNode',1, ...
% %          'LabelNode', 11, 'MarkerSize',4);
% %     PlotSBFEMesh(d{id}.coord, d{id}.sdConn, opt); 
% %     figure
end
%ltx resolution of coordinates (round to the first digit in 1, 2 or 5). 
res = 0.1*minLeng; %ltx the factor is arbitrarily chosen as 0.1 \label{createMultiDomainMesh_eps}
expnt = floor(log10(res)); 
firstDigit = floor(res/10^expnt);
roundTo = [1 2 2 5 5 5 5 10 10];
res = 10^expnt*roundTo(firstDigit); %ltx round the resolution
%ltx Combining the domains one by one iteratively \label{createMultiDomainMeshCombine}
da = d{1}; %ltx including the first domain
for id = 2:nDomain; %ltx loop over the remaining domains
    [ da.coord, da.sdConn, da.sdSC ] = combineTwoMeshes( ...
        {d{id}, da}, dDefs.onLine{id-1}, res);
    [ da.meshEdge, da.sdEdge] = meshConnectivity( da.sdConn );
end

coord = da.coord; sdConn = da.sdConn; sdSC = da.sdSC; %ltx output

end