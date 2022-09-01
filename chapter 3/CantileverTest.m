function [Res]=CantileverTest()
%---------------------------------------------------------------------input
Lx = 48;   Ly = 12;   h  = 1;   E  = 30e03; nu = 0.25;
Res = zeros(5,2);
for i = 1:size(Res,1)
sx      = 2^(i+1);           sy      = 2^i;   
%---------------------------------------------------------------------nodes
x       = 0:Lx/sx:Lx;        y       = 0:Ly/sy:Ly;    
[X,Y]   = meshgrid(x,y);
num     = size(X(:),1);      nID     = (1:num)';
Nodes   = [nID,X(:),Y(:),0*X(:)];
%------------------------------------------------------------------elements
Trigs   = delaunay(round(X(:),3),round(Y(:),3));
%----------------------------------------------------------------------dofs
Xcon    = nID(X==0);
rDofs   = [Xcon,0*Xcon+1; Xcon,0*Xcon+2;Xcon,0*Xcon+6];
AllDof  = reshape((1:6*num),6,num);
ActDof  = AllDof([1,2,6],:);
NonDof  = 6*(rDofs(:,1)-1)+rDofs(:,2);
dofs    = setdiff(ActDof(:),NonDof);
%----------------------------------------------------------------------load
vw      = 6*(nID((X==Lx)&(Y==Ly/2))-1)+2;
vw2      = 6*(nID(abs(X-Lx)<1e-5)-1)+2;
P       = zeros(6*num,1);    P(vw2,1)=40/size(vw2(:),1);
%-----------------------------------------------------------------------FEM
U       = 0*P;
K       = KmembrMaker(Nodes,Trigs,h,E,nu);
U(dofs) = K(dofs,dofs)\P(dofs);
Res(i,2)  = U(vw)/(40*Lx^3/(3*E*(h*Ly^3/12)));
Res(i,1) = sx;
end
figure
plot(Res(:,1),Res(:,2))
xlabel('mesh size') 
ylabel('displacement sensitivity') 
end

function [K] = KmembrMaker(Nodes,Trig,h, E, nu)
K = zeros(6*size(Nodes,1));
for i = 1:size(Trig,1)
    Cor  = Nodes(Trig(i,:),2:4);
    Te   = TeMatrix(Cor);
    Tmat = kron(eye(3),[Te(1:2,:),zeros(2,3);zeros(1,3),Te(3,:)]);
    locC = (Cor-repmat(Cor(1,:),3,1))*Te';
    KeL  = KeLCS_TrigMembraneDrill(locC(:,1:2), h, E, nu);
    Ke   = Tmat'*KeL*Tmat;
    pos = [6*(Trig(i,1)-1)+(1:6) 6*(Trig(i,2)-1)+(1:6) 6*(Trig(i,3)-1)+(1:6)];
    K(pos,pos) = K(pos,pos)+Ke;
end
end

function [Te]=TeMatrix(Coords)
%-----------------------------------------------------------------------e_z
NumerZ  = cross(Coords(2,:)-Coords(1,:),Coords(3,:)-Coords(1,:));
eZ      = NumerZ/norm(NumerZ);
%-----------------------------------------------------------------------e_x
NumerX  = 0.5*(Coords(2,:)+Coords(3,:))- Coords(1,:);
eX      = NumerX/norm(NumerX);
%-----------------------------------------------------------------------e_y
NumerY  = cross(eZ,eX);
eY      = NumerY/norm(NumerY);
Te = [ eX; eY; eZ];
end

