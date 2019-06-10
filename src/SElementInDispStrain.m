function [nodexy, dsp, strnNode, GPxy, strnEle] = ...
    SElementInDispStrain(xi, sdSln, sdStrnMode, sdIntgConst)
%Displacements and strains at specified radial coordinate
% Inputs:
%   xi            - radial coordinate
%   sdSln         - solutions for S-element 
%   sdStrnMode    - strain modes of S-element
%   IntgConst     - integration constants
%
% Outputs (All valus are on the scaled boundary, i.e. coodinate
%          line, at specified xi):
%   nodexy(i,:)   - coordinates of node i
%   dsp(i,:)      - nodal displacement funcitons of node i 
%   strnNode(:,i) - strains on the radial line of node i
%   GPxy(ie,:)    - coordinates of middle point of element ie
%   strnEle(:,ie) - strains at middle point of element ie

%cbgnltx
% Transform local coordinates (origin at scaling centre) 
% at scaled boundary to global coordinates.\\
% \texttt{GPxy(ie,:)} -  coordinates of Gauss point 
% (middle of 2-node element) of  element \texttt{ie} after scaling.\\
% \texttt{nodexy(i,:)} - coordinates of node \texttt{i} after scaling.
%cendltx
GPxy   = bsxfun(@plus, xi*sdStrnMode.xy, sdSln.sc); 
nodexy = bsxfun(@plus, xi*sdSln.xy, sdSln.sc);

if xi >1.d-16 %ltx outside of a tiny region around the scaling centre
    fxi = (xi.^sdSln.d).*sdIntgConst; %ltx $\xi^{\left\langle \lambda_{b}\right\rangle }\{c\}$
    dsp = sdSln.v*fxi; %ltx Eq.~\eqref{eq:EigenMethod-matrix-u}
    strnEle = sdStrnMode.value(:,1:end-2) ...
               *fxi(1:end-2)/xi; %ltx Eq.~\eqref{eq:EignMethod-strain-matrix}
else %ltx at scaling centre
    dsp = sdSln.v(:,end-1:end)* ...
          sdIntgConst(end-1:end); %ltx Eq.~\eqref{eq:EignMethod-disp-scalingCentre}
    if(min(real(sdSln.d(1:end-2)))>0.999)
        strnEle = sdStrnMode.value(:,end-5:end-2)...
            *sdIntgConst(end-5:end-2); %ltx Eq.~\eqref{eq:EignMethod-strain-scalingCentre}
    else %ltx stress singularity at scaling centre
        strnEle = NaN(length(sdStrnMode.value),1); 
    end
end
%ltx remove possible tiny imaginary part due to numerical error 
dsp = real(dsp);
strnEle = real(strnEle);

%ltx \texttt{strnEle(1:3,ie)} is the strains at centre of element \texttt{ie} after reshaping.
strnEle = reshape(strnEle, 3, []);

%ltx nodal stresses by averaging element stresses
nNode = length(sdSln.node); %ltx number of nodes
strnNode = zeros(3,nNode); %ltx initialisation
%ltx counters of number of elements connected to a node
count = zeros(1,nNode);
n = sdSln.conn(:,1); % vector of first node of elements
%ltx add element stresses to first node of elements
strnNode(:,n) = strnNode(:,n) + strnEle; 
%ltx increment counters
count(n) = count(n) + 1;

n = sdSln.conn(:,2); % vector of second node of elements
%ltx add element stresses to second node of elements
strnNode(:,n) = strnNode(:,n) + strnEle;
%ltx increment counters
count(n) = count(n) + 1;
strnNode = bsxfun(@rdivide, strnNode, count); %ltx averaging

end