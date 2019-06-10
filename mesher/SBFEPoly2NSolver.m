function varargout = SBFEPoly2NSolver(probdef, ...
                               coord, sdConn, sdSC, opt)
%cbgnltx 
%\myCodeTitle{% 
%Scaled boundary finite element analysis using polygon
% S-elements discretised with 2-node line elements at boundary}
%{\bf Inputs:}
%\begin{myDescription}
%\item[probdef]  handle of function defining the problem
%\item[coord(i,:)]   coordinates of node \texttt{i}
%\item[sdConn\{isd,:\}(ie,:)]  S-element conncetivity. The nodes of
%                         line element \texttt{ie} in S-element \texttt{isd}.
% \item[sdSC(isd,:)]  coordinates of scaling centre of S-element
% \texttt{isd}
% \item[opt] options of analysis. The fields are:\\
% \begin{myTabular}{lll}
%   \texttt{type} & - &  type of analysis, =\texttt{'STATICS'}, 
%                         \texttt{'MODAL'} or \texttt{'TIME'}\\
%   \texttt{modalPara}& - & number of modes \\
%   \texttt{TIMEPara} & - & an array with [ number of time steps, size of time step]
% \end{myTabular}
%\end{myDescription}
%{\bf Outputs:}\par
%\begin{enumerate}
%	\item \emph{For static analysis}\\
%\begin{myTabular}{lll}
% \texttt{U} & - &  nodal displacements\\
% \texttt{sdSln} & - &  solutions for S-element
% \end{myTabular}
%
%\item \emph{For modal analysis}\\
%\begin{myTabular}{lll}
% \texttt{freqAngular}  & - &  angular frequencies\\
% \texttt{modeShape} &  - &   modal shapes
% \end{myTabular}
%
%\item \emph{For response history analysis}\\
%\begin{myTabular}{lll}
% \texttt{tn} & - & time stations \\
% \texttt{dsp} & - & nodal displacements \\
% \texttt{vel} & - & nodal velocities \\
% \texttt{accl} & - & nodal accelerations
% \end{myTabular}
%\end{enumerate}
%cendltx

%ltx {\bf Material}: elasticity matrix and mass density\label{SBFEPoly2NMaterial}
mat = probdef('MAT'); %ltx get material constants

%ltx {\bf S-element solution and global stiffness and mass matrix}\label{SBFEPoly2NSubDm}
[sdSln, K, M] = SBFEMAssembly(coord, sdConn, sdSC, mat);

if nargin < 5
    opt.type = 'STATICS';
end

switch opt.type
    case 'STATICS'; %ltx \bf{Static analysis} \label{SBFEPoly2NStatics}
        %ltx \emph{Boundary conditions} 
        %ltx \verb!outArg{1}!: displacement constrains. \verb!outArg{2}!: nodal forces  
        outArg = probdef('BCond', {coord, sdConn});
        %ltx solution of nodal displacements and reactions \label{SBFEPoly2NSolution}
        [U, ReactFrc] = SolverStatics(K, outArg{1}, outArg{2});
        
        %ltx \emph{Post-processing} \label{SBFEPoly2NPostProcessing}
        %ltx plot deformed mesh
        figure
        opt = struct('MagnFct', 0.1, 'Undeformed','--k');
        PlotDeformedMesh(U, coord, sdConn, opt);
        title('DEFORMED MESH');
        %ltx output displacements at selected DOFs
        output = probdef('Output', {coord, sdConn});
        if ~isempty(output)
            disp('   DOF          Displacement');
            fprintf('%6d  %25.15e\n',[output.DispDOF ...
                                      U(output.DispDOF)]');
        end
        %ltx compute error norm
        Uex = probdef('EXACT',coord);
        if ~isempty(Uex)
            %ltx reshape the exact solution to a column vector
            Uex = reshape(Uex(:,1:2)',[],1);
            err = Uex-U; %ltx displacement error
            eNorm = norm(err)/norm(Uex);
            disp([' Displacement error norm = ', num2str(eNorm)]);
        end
        %ltx output
        varargout = {U, sdSln};
        
    case 'MODAL' %ltx {\bf Modal analysis} \label{SBFEPoly2NModal}
        NMode = opt.modalPara; %ltx number of modes
        outArg = probdef('BCondModal', {coord, sdConn});
        [freqAngular, modeShape] = SolverMode(NMode, K, M, ...
               outArg{1}); %ltx natural frequencies and modal shapes
        %ltx output natrual frequency and period  
        fprintf('%6s %17s %17s\n', 'Mode', 'Frequency', 'Period')
        for imode = 1:NMode; %ltx mode number
            fprintf('%6d   %15.7d   %15.7d\n', imode, ...
                freqAngular(imode), 2*pi./freqAngular(imode));
            figure %ltx plot modal shape
            opt = struct('MagnFct', 0.1, 'Undeformed','--k');
            PlotDeformedMesh(modeShape(:,imode),coord,sdConn,opt);
            title(['MODE SHAPE: ',int2str(imode)]);
        end
        varargout = {freqAngular, modeShape};
        
    case 'TIME' %ltx {\bf Response history analysis} \label{SBFEPoly2NTime}
        TIMEPara = opt.TIMEPara;
        outArg = probdef('BCondTime', {coord, sdConn});
        [tn, dsp, vel, accl] = SolverNewmark(TIMEPara, ...
                           outArg{3}, K, M, outArg{1}, outArg{2});
        %ltx plot displacements at selected DOFs
        output = probdef('OutputTime', {coord, sdConn});
        if ~isempty(output)
            figure
            plot(tn,dsp(output.DispDOF,:),'-k')
            xlabel('Time (sec)');
            ylabel('Displacment (m)');
% %             title(['DISPLACEMENT AT DOFs: ',num2str(output.DispDOF')]);
            title('DISPLACEMENT RESPONSE');
        end
        varargout = {tn, dsp, vel, accl};
    otherwise
        
end