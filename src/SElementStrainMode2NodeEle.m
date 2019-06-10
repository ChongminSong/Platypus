function [ sdStrnMode ] = SElementStrainMode2NodeEle( sdSln )
%Interal displacements and strain modes of S-elements
%
% Inputs:
%   sdSln       - solutions for S-elements 
% Outputs:
%   sdStrnMode  - strain modes of S-elements

Nsd = length(sdSln); %ltx number of S-elements
sdStrnMode = cell(Nsd,1); %ltx initialisation of output argument
for isd = 1:Nsd %ltx loop over S-elements
    
    %ltx number of DOFs of boundary nodes
    nd = 2*length(sdSln{isd}.node); 
    v = sdSln{isd}.v; %ltx displacement modes
   
    %ltx {\bf Strain modes}. The last 2 columns equal 0 (rigid body motions)
    %ltx number of DOFs excluding 2 translational rigid-body motions
    nd2 = nd-2; 
    d = sdSln{isd}.d(1:nd2); %ltx eigenvalues
    %cbgnltx
    % See Eq.~\eqref{eq:EignMethod-strainModes-matrix}. 
    % $[\Phi_{b}^{(\mathrm{u})}]\left\langle \lambda_{b}\right\rangle$ is
    % computed. The element values are extracted based on element
    % connectivity in Line~\ref{SElementStrainMode2NodeEleStrnMode}
    %cendltx
    vb = v(:,1:nd2) .* (d(1:nd2,ones(1,nd)))'; 
    
    n1 = sdSln{isd}.conn(:,1); %ltx first node of all 2-node elements
    n2 = sdSln{isd}.conn(:,2); %ltx second node of all 2-node elements
    %ltx \texttt{LDof(i,:)}: Local DOFs of nodes of element \texttt{i} in an S-element
    LDof = [ n1+n1-1 n1+n1 n2+n2-1 n2+n2 ];
    xy = sdSln{isd}.xy; %ltx nodal coordinates with origin at scaling centre
    %ltx \texttt{dxy(i,:)}: $[\Delta_{x},\ \Delta_{y}]$ of $i$-th element, Eq.~\eqref{eq:2D-2NodeElement-dxy}
    dxy = xy(n2,:) - xy(n1,:);
    %ltx \texttt{mxy(i,:)}:  $[\bar{x},\ \bar{y}]$ of $i$-th element, Eq.~\eqref{eq:2D-2NodeElement-avgxy}
    mxy = (xy(n2,:) + xy(n1,:))/2;
    %ltx \texttt{a(i)}:  $2|J_b|$ of $i$-th element, Eq.~\eqref{eq-Jdet-1}
    a = mxy(:,1).*dxy(:,2) - mxy(:,2).*dxy(:,1);
    
    ne = length(n1); %ltx numer of line elements
    mode = zeros(3*ne, nd); %ltx initialising strain modes
    for ie = 1:ne %ltx loop over elements at boundary
        C1 = 0.5*[ dxy(ie,2) 0; 0 -dxy(ie,1); ...
            -dxy(ie,1) dxy(ie,2)];  %ltx Eq.~\eqref{eq:2D-2NodeElement-C1}
        C2 = [-mxy(ie,2) 0; 0 mxy(ie,1); ...
            mxy(ie,1) -mxy(ie,2)]; %ltx Eq.~\eqref{eq:2D-2NodeElement-C2}
        B1 = 1/a(ie)*[ C1 C1]; %ltx Eq.~\eqref{eq:2NodeEle-GP-B1}
        B2 = 1/a(ie)*[-C2 C2]; %ltx Eq.~\eqref{eq:2NodeEle-GP-B2}
        mode(3*(ie-1)+1:3*ie,1:nd2) = B1*vb(LDof(ie,:),:) ...
            + B2*v(LDof(ie,:),1:nd2);  %ltx strain modes, Eq~\eqref{eq:EignMethod-strainModes-matrix} \label{SElementStrainMode2NodeEleStrnMode}
    end
    
    %cbgnltx
    %\increaseIndent Store the ouput in cell array \texttt{sdPstP}.
    % The number of S-element \texttt{isd} is the index of the array.\\
    % \texttt{GPxy(ie,:)}: the coordinates of the Gauss Point
    % of element \texttt{ie} (middle point of 2-node element).\\
    % \texttt{strnMode(:,ie)}: the strain modes at the Gauss Point
    % of element \texttt{ie}.
    %cendltx
    sdStrnMode{isd}= struct('xy',mxy, 'value',mode);
    
end