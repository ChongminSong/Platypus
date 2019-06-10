function varargout = SBFEPolyStress(probdef, U, sdSln, ...
                               sdConn, sdSC) 
%compute stress of each node and scaling centre of each polygon subdomain
%
% Inputs:
%   U                 - nodal displacements
%   sdSln{isd,:}      - solutions for subdomain isd
%   sdConn{isd,:}     - conncetivity of subdomain isd
%   sdSC(isd,:)       - coordinates of scaling centre of subdomain isd
%
% Outputs:
%   SCStrs(i,:)       - stress of scaling center of polygon subdomain i
%   NodeStrs          - nodal stress                             

%ltx {\bf Material}: elascity matrix and mass density\label{SBFEPoly2NMaterial}
       mat = probdef('MAT'); %ltx get material constants
%%...compute strain mode and integration constants of each subdomain
       sdStrnMode  = SElementStrainMode2NodeEle(sdSln);
       sdIntgConst = SElementIntgConst(U, sdSln);
       
       n=size(sdSC,1); 
       %%initialisation of matrix of scaling centre stress 
       SCStrs=zeros(n,3);
       %%initialisation of nodal stress, NStrs(:,4) store
       %%counts for averaging  
       NStrs=zeros(size(U,1)/2,4);
      
       for i=1:n  
        %%...compute Scaling Centre strain, xi==0  
        [~, ~, strnNode0, ~, strnEle0]=SElementInDispStrain(0,sdSln{i},  ...
        sdStrnMode{i}, sdIntgConst{i});
        %%...compute nodal strain on the boundary, xi==1
        [~, ~, strnNode1, ~, strnEle1]=SElementInDispStrain(1,sdSln{i},  ...
        sdStrnMode{i}, sdIntgConst{i});
        %%...compute Scaling Centre stress
        strsEle0=mat.D*strnEle0(:,1);
        SCStrs(i,1:3)=strsEle0';
        %%...compute nodal stress on the boundary 
        strsNode1=mat.D*strnNode1;
        %%...left column for counts
        strsNode1=[strsNode1' ones(size(strsNode1,2),1)];
        NStrs(sdConn{i}(:,1),:)=NStrs(sdConn{i}(:,1),:)+strsNode1;
       end
       
       %%...average nodal stress
          NodeStrs=NStrs(:,1:3)./repmat(NStrs(:,4),[1 3]);
        
       %%...compute scaling centre stress error norm
            S=reshape(SCStrs(:,1:3), [], 1);
            SCStrsEx = probdef('EXACT',sdSC); 
         if ~isempty(SCStrsEx)
            SCStrsEx=reshape(SCStrsEx(:,3:5), [], 1);
            Serr = SCStrsEx-S; %ltx Scaling Centre Stress error
            SNorm = norm(Serr)/norm(SCStrsEx);
            disp([' Scaling Centre Stress error norm = ', num2str(SNorm)]);
         end
        
        varargout = {SCStrs, NodeStrs};           
end
        
