function [ e0, e1, e2, m0 ] = EleCoeff2NodeEle( xy, mat)
%Coefficient matrices of a 2-node line element
%
%  Inputs:
%   xy(i,:)    - coordinates of node i (orgin at scaling centre).
%                The nodes are numbered locally within
%                each line element
%   mat        - material constants
%     mat.D    - elasticity matrix
%     mat.den  - mass density
%
%  Outputs:
%    e0, e1, e2, m0  - element coefficient matrices

dxy = xy(2,:)-xy(1,:); %ltx $[\Delta_{x},\ \Delta_{y}]$, Eq.~\eqref{eq:2D-2NodeElement-dxy} 
mxy = sum(xy)/2; %ltx $[\bar{x},\ \bar{y}]$, Eq.~\eqref{eq:2D-2NodeElement-avgxy}
a = xy(1,1)*xy(2,2)-xy(2,1)*xy(1,2); %ltx $a=2|J_b|$, Eq.~\eqref{eq-Jdet}
if a < 1.d-10
    disp('negative area (EleCoeff2NodeEle)');
    pause
end
C1 = 0.5*[ dxy(2) 0; 0 -dxy(1); -dxy(1) dxy(2)];  %ltx Eq.~\eqref{eq:2D-2NodeElement-C1}
C2 = [-mxy(2) 0; 0 mxy(1); mxy(1) -mxy(2)]; %ltx Eq.~\eqref{eq:2D-2NodeElement-C2}
Q0 = 1/a*(C1'*mat.D*C1); %ltx Eq.~\eqref{eq:2D-2NodeElement-Q0}
Q1 = 1/a*(C2'*mat.D*C1); %ltx Eq.~\eqref{eq:2D-2NodeElement-Q1}
Q2 = 1/a*(C2'*mat.D*C2); %ltx Eq.~\eqref{eq:2D-2NodeElement-Q2} 
e0 =  2/3*[2*Q0 Q0; Q0 2*Q0]; %ltx  Eq.~\eqref{eq:2D-2NodeElement-E0}
e1 = -1/3*[ Q0 -Q0; -Q0  Q0] + [-Q1 -Q1;  Q1 Q1]; %ltx  Eq.~\eqref{eq:2D-2NodeElement-E1}
e2 =  1/3*[ Q0 -Q0; -Q0  Q0] + [ Q2 -Q2; -Q2 Q2]; %ltx  Eq.~\eqref{eq:2D-2NodeElement-E2}

%ltx mass coefficent matrix, Eq.~\eqref{eq:2D-2NodeElem-M0} \label{Code_ElemCoeff_mass}
m0 = a*mat.den/6*[ 2 0  1 0; 0 2  0 1; 1 0  2 0; 0 1 0 2 ];

end