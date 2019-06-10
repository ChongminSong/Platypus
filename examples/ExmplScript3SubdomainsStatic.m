% input problem definition
Exmpl3Subdomains

% solution of subdomains and global stiffness and mass matrices
[sdSln, K, M] = SBFEMAssemblage(coord, sdConn, sdSC, mat);

% static solution
dofs = (1:numel(coord))'; % vector of DOFs for output 
SBFEM2DStaticsScript
