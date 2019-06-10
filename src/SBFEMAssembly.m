function [sdSln, K, M] = SBFEMAssembly(coord, sdConn, sdSC, mat)
%Assembly of global stiffness and mass matrices 
%
%Inputs:
%  coord(i,:)   - coordinates of node i
%  sdConn{isd,:}(ie,:)  - S-element conncetivity. The nodes of 
%                         line element ie in S-element isd. 
%  sdSC(isd,:)  - coordinates of scaling centre of S-element isd
%  mat         - material constants
%     mat.D     - elasticity matrix
%     mat.den   - mass density
%
%Outputs:
%  sdSln        - solutions for S-element 
%  K            - global stiffness matrix 
%  M            - global mass matrix 

%ltx {\bf Solution of subdomains}
Nsd = length(sdConn); %ltx number of subdomains
sdSln = cell(Nsd,1); %ltx store solutions for subdomains
for isd = 1:Nsd %ltx loop over subdomains
    %cbgnltx
    % \increaseIndent
    % \texttt{sdNode} contains global nodal numbers of the nodes in  
    % an S-element. Vector \texttt{ic} maps the global connectivity to 
    % the local connectivity of the S-element
    %cendltx
    [sdNode, ~, ic] = unique(sdConn{isd}(:)); %ltx  remove duplicats \label{Code_LocalConn}
    sdNode = reshape(sdNode,1,[]); %ltx nodes stored in one row
    xy = coord(sdNode,:); % nodal coordinates
    %ltx transform coordinate origin to scaling centre\label{Code_LocalCoord}
    xy = bsxfun(@minus, xy, sdSC(isd,:)); 
    %ltx line element connectivity in local nodal numbers of an S-element
    LConn = reshape(ic, size(sdConn{isd}) );
    %ltx compute S-element coefficient matrices
    [ E0, E1, E2, M0 ] = SElementCoeffMtx(xy, LConn, mat);
    %ltx compute solution for S-element
    [ K, d, v, M] = SElementSlnEigenMethod(E0, E1, E2, M0);
    
    %ltx store S-element data and solution \label{Code_StoreSubdomainSln}
     sdSln{isd}= struct('xy',xy, 'sc',sdSC(isd,:),  ...
         'conn',LConn, 'node',sdNode, ...
         'K',K, 'M',M, 'd',d, 'v', v);
end

%ltx {\bf Assembly} \label{Code_assemblage}
%ltx sum of entries of stiffness matrices of all subdomains
ncoe = sum(cellfun(@(x) numel(x.K), sdSln));
%ltx initialising non-zero entries in global stiffness and mass matrix
K = zeros(ncoe,1); M = zeros(ncoe,1); 
%ltx rows and columns of non-zero entries in global stiffness matrix
Ki = K; Kj = K; 
StartInd = 1; %ltx starting position of an S-element stiffness matrix
for isd = 1:length(sdSln); %ltx loop over subdomains
    %ltx global DOFs of nodes in an S-element \label{Code_subdomainDOF}
    dof = reshape([2*sdSln{isd}.node-1;2*sdSln{isd}.node], [], 1);
    Ndof = length(dof); %ltx number of DOFs of an S-element
    %ltx rRow and column numbers of stiffness coefficients of an S-element
    sdI = repmat(dof,1, Ndof); sdJ = sdI'; 
    
    %ltx store stiffness, row and column indices
    EndInd = StartInd + Ndof^2 - 1; %ltx ending position 
    K (StartInd:EndInd) = sdSln{isd}.K(:);
    M (StartInd:EndInd) = sdSln{isd}.M(:);
    Ki(StartInd:EndInd) = sdI(:);
    Kj(StartInd:EndInd) = sdJ(:);
    
    StartInd = EndInd + 1;  %ltx increment the starting position
end
%ltx form global stiffness matrix
K = sparse(Ki,Kj,K); %ltx form global stiffness matrix in sparse storage
K = (K+K')/2; %ltx ensure symmetry
%ltx form global mass matrix
M = sparse(Ki,Kj,M); %ltx Form global mass matrix in sparse storage
M = (M+M')/2; %ltx ensure symmetry

end