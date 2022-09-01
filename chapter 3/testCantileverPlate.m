power = (-6:-6)';
nTest = size(power,1);
resSA = cell(nTest,1);
resFD = cell(nTest);
for i=1:nTest
    step    = 10^power(i);
    history = comparMethods(step);
    resSA{i} = history.FD;
    resFD{i} = history.Fg;
end
resSA = cell2mat(resSA);
resFD = cell2mat(resFD);
function [history] = comparMethods(step)
%% VARIABLES
% geometry
Lx      = 1;            % length of plate
Ly      = 1;            % width of plate
t       = 0.1;          % thickness of plate
% material
E       = 1.09e6;
nu      = 0.3;

for ind = 0:4
    f       = 2^ind;
    fx      = f;
    fy      = 2*f;
    %% DEFINITION OF THE NODES
    num     = (fx+1)*(fy+1); 
    nID     = (1:num)';
    [X,Y,Z] = meshgrid((0:fx)*Lx/fx,(0:fy)*Ly/fy,0);
    Nodes   = [nID,X(:),Y(:),Z(:)];                                                  % difinition of the nodes
    %% DEFINITION OF THE ELEMENTS
    Trigs   = delaunay(X,Y);
    %% DEFINITION OF THE DOFs
    free    = nID(X(:)~=0);                      
    dofTZ   = 3*(free-1)+1;
    dofRX   = 3*(free-1)+2;
    dofRY   = 3*(free-1)+3;
    dofs    = [dofTZ; dofRX; dofRY];
    %% DEFINITION OF THE LOADS
    IDs     = nID(X(:)==Lx);
    load    = conv2(ones(fy,1)/fy,[0.5;0.5]);
    nDOF    = 3*num;
    P       = sparse(3*(IDs-1)+1,0*IDs+1,load,nDOF,1);
    %% DEFINITION OF THE STIFFNESS MATRIX
    duAN    = Analytical(Nodes,Trigs,t,E,nu,P,dofs);
    duFD    = FiniteDiff(Nodes,Trigs,t,E,nu,P,dofs,step);
    duFg    = FiniteDiffGlob(Nodes,Trigs,t,E,nu,P,dofs,step,fx,fy);
    history.subd(:,ind+1) = f;
    indEnd  = 3*((fx)*(fy+1)+1:(fx+1)*(fy+1))-2;
    history.FD(:,ind+1) = mean(sum(duFD(indEnd,:),2));
    history.Fg(:,ind+1) = mean(sum(duFg(indEnd,:),2));
    history.an(:,ind+1) = mean(sum(duAN(indEnd,:),2));

   
end


figure
plot(history.subd,history.FD,history.subd,history.Fg,':',history.subd,history.an,'--')
xlabel('mesh size') 
ylabel('Normalized tip deflection') 

end


 %________________________________________________Fin_diff_Stiffness_Matrix
function [du] = FiniteDiff(Nodes,Trigs,h,E,nu,P,dofs,step)
     U  = 0*P;    
     K  = GenerateK(Nodes,Trigs,h,E,nu);
U(dofs) = K(dofs,dofs)\P(dofs);
ids = Nodes(Nodes(:,2)==max(Nodes(:,2)),1);
du  = zeros(size(P,1),size(ids,1));
for i   = 1:size(ids,1)
    dKdx        = GenerateFinDiffK(Nodes,Trigs,h,E,nu,ids(i),step);
    du(dofs,i)  = K(dofs,dofs)\(-dKdx(dofs,dofs)*U(dofs));
end
end
function [dKdx] = GenerateFinDiffK(Nodes,Trigs,h,E,nu,i,step)
UseTrig  = Trigs((Trigs(:,1)==i)|(Trigs(:,2)==i)|(Trigs(:,3)==i),:);
dKdx     = FinDiffKAll(Nodes,UseTrig,h,E,nu,i,step);
end
function [dKdx] = FinDiffKAll(Nodes,Trigs,h,E,nu,i,step)
nDOF = 3*size(Nodes,1);      nTRI = size(Trigs,1);
  Mi = zeros(81*nTRI,1);      Mj = Mi;  Mz = Mi;
for j=1:nTRI
      Coor  = Nodes(Trigs(j,:),2:4);
         d  = sum((Trigs(j,:)==i)*(1:3)',1);
     dKedx  = FinDiffKeGCS_TrigPlateDKT_Coor(Coor,h,E,nu,d,step);
       pos  = Addkron(3*(Trigs(j,:)-1),(1:3));
       ind  = 81*(j-1)+(1:81);
    Mi(ind) = kron(ones(1,9),pos);
    Mj(ind) = kron(pos,ones(1,9));
    Mz(ind) = dKedx(:);
end
dKdx = sparse(Mi,Mj,Mz,nDOF,nDOF);
end
function [dKedx]=FinDiffKeGCS_TrigPlateDKT_Coor(Coords,h,E,nu,d,step)
Te      = kron(eye(3),[1 0 0; 0 0 -1; 0 1 0]);
     CL  = Coords;
     CA  = Coords; CA(d)=CA(d)+step;
     CLP = CA;
    KeLP = KeLCS_TrigPlateDKT(CLP(:,1:2), h, E, nu);
    KeL  = KeLCS_TrigPlateDKT(CL(:,1:2), h, E, nu);
   dKeL  = (KeLP-KeL)/step;
  dKedx  = Te'*dKeL*Te;
end





function [du] = FiniteDiffGlob(Nodes,Trigs,h,E,nu,P,dofs,step,fx,fy)
     U  = 0*P;    
     K  = GenerateK(Nodes,Trigs,h,E,nu);
U(dofs) = K(dofs,dofs)\P(dofs);
ids = (fx*(fy+1)+1:(fx+1)*(fy+1));
du  = zeros(size(P,1),1);
    NodP             = Nodes; 
    NodP(ids,2)   = NodP(ids,2)+step;
    Kper             = GenerateK(NodP,Trigs,h,E,nu);
    dKdx             = (Kper-K)/step;
    du(dofs,:)       = K(dofs,dofs)\(-dKdx(dofs,dofs)*U(dofs));
end 




function [du] = Analytical(Nodes,Trigs,h,E,nu,P,dofs)
     U  = 0*P;    
     K  = GenerateK(Nodes,Trigs,h,E,nu);
U(dofs) = K(dofs,dofs)\P(dofs);
ids = Nodes(Nodes(:,2)==max(Nodes(:,2)),1);
du  = zeros(size(P,1),size(ids,1));
for i   = 1:size(ids,1)
    dKdx        = GenerateDiffK(Nodes,Trigs,h,E,nu,ids(i));
    du(dofs,i)  = K(dofs,dofs)\(-dKdx(dofs,dofs)*U(dofs));
end
end
%__________________________________________________________Stiffness_Matrix
function [K] = GenerateK(Nodes,Trigs,t,E,nu)
nDOF = 3*size(Nodes,1);      nTRI = size(Trigs,1);
  Mi = zeros(81*nTRI,1);      Mj = Mi;  Mv = Mi;
for i = 1:nTRI
    Cor     = Nodes(Trigs(i,:),2:4);
    Te      = kron(eye(3),[1 0 0; 0 0 -1; 0 1 0]);
    locC    = (Cor-repmat(Cor(1,:),3,1));
    KPl     = KeLCS_TrigPlateDKT(locC(:,1:2),t,E,nu);
    Ke      = Te'*KPl*Te;
    pos     = Addkron(3*(Trigs(i,:)-1),(1:3));
    ind     = 81*(i-1)+(1:81);
    Mi(ind) = kron(ones(1,9),pos);
    Mj(ind) = kron(pos,ones(1,9));
    Mv(ind) = Ke(:);
end
K = sparse(Mi,Mj,Mv,nDOF,nDOF);
end



%_____________________________________________________Diff_Stiffness_Matrix
function [dKdx] = GenerateDiffK(Nodes,Trigs,h,E,nu,i)
UseTrig  = Trigs((Trigs(:,1)==i)|(Trigs(:,2)==i)|(Trigs(:,3)==i),:);
dKdx     = diffKAll(Nodes,UseTrig,h,E,nu,i);
end
function [dKdz] = diffKAll(Nodes,Trigs,h,E,nu,i)
nDOF = 3*size(Nodes,1);      nTRI = size(Trigs,1);
  Mi = zeros(81*nTRI,1);      Mj = Mi;  Mz = Mi;
for j=1:nTRI
      Coor  = Nodes(Trigs(j,:),2:4);
         d  = sum((Trigs(j,:)==i)*(1:3)',1);
     dKedz  = DiffKeGCS_TrigPlateDKT_Coor(Coor,h,E,nu,d);
       pos  = Addkron(3*(Trigs(j,:)-1),(1:3));
       ind  = 81*(j-1)+(1:81);
    Mi(ind) = kron(ones(1,9),pos);
    Mj(ind) = kron(pos,ones(1,9));
    Mz(ind) = dKedz(:);
end
dKdz = sparse(Mi,Mj,Mz,nDOF,nDOF);
end


% PLATE STIFFNESS MATRIX %-------------------------------------------------
function [dKedx]=DiffKeGCS_TrigPlateDKT_Coor(Coords,h,E,nu,d)
Te      = kron(eye(3),[1 0 0; 0 0 -1; 0 1 0]);
dTe     = 0*Te;
dCoords  = d==reshape((1:9),3,3);
     CL  =  Coords-repmat( Coords(1,:),3,1);
    dCL  = dCoords-repmat(dCoords(1,:),3,1);
           [dKdx1,dKdy1,dKdx2,dKdy2,dKdx3,dKdy3] = ...
           DiffKeLCS_TrigPlateDKT_Coor(CL(:,1:2), h, E, nu);
    KeL  = KeLCS_TrigPlateDKT(CL(:,1:2), h, E, nu);
   dKeL  = dKdx1*dCL(1,1)+dKdx2*dCL(2,1)+dKdx3*dCL(3,1)+...
           dKdy1*dCL(1,2)+dKdy2*dCL(2,2)+dKdy3*dCL(3,2);
  dKedx  = dTe'*KeL*Te+Te'*dKeL*Te+Te'*KeL*dTe;
end
