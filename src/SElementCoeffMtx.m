function [ E0, E1, E2, M0 ] = SElementCoeffMtx(xy, conn, mat)
%Coefficient matrix of an S-element
%
% Inputs:
%   xy(i,:)    - coordinates of node i (orgin at scaling centre) 
%                The nodes are numbered locally within 
%                an S-element starting from 1
%   conn(ie,:) - local connectivity matrix of line element ie 
%                in the local nodal numbers of an S-element
%   mat        - material constants
%     mat.D    - elasticity matrix
%     mat.den  - mass density
%
% Outputs:
%   E0, E1, E2  - coefficient matrices of S-element

nd = 2*size(xy,1); % number of DOFs at boundary (2 DOFs per node)
%ltx initialising variables
E0 = zeros(nd, nd);
E1 = zeros(nd, nd);
E2 = zeros(nd, nd);
M0 = zeros(nd, nd);
for ie = 1:size(conn,1) %ltx loop over elements at boundary
    xyEle = xy(conn(ie,:),:); %ltx nodal coordinates of an element \label{Code_subdomainCoeff_nodalcoord} 
    %ltx get element coefficient matrices of an element 
    [ ee0, ee1, ee2, em0 ] = EleCoeff2NodeEle(xyEle, mat);
    %ltx local DOFs (in S-element) of an element\refstepcounter{lineno}\label{Code_subdomainCoeff_conn}\addtocounter{lineno}{-1}
    d = reshape([2*conn(ie,:)-1; 2*conn(ie,:)], 1, []);
    %ltx assemble coefficient matrices of S-element\refstepcounter{lineno}\label{Code_subdomainCoeff_E0}\addtocounter{lineno}{-1}
    E0(d,d) = E0(d,d) + ee0; 
    E1(d,d) = E1(d,d) + ee1;
    E2(d,d) = E2(d,d) + ee2;
    M0(d,d) = M0(d,d) + em0;
end

end