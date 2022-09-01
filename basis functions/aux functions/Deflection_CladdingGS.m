function [fC] = Deflection_CladdingGS(Nodes,Trigs,Pdis,varargin)
if nargin < 4
	C  = CompCladding(Nodes,Trigs,Pdis);
    dNd = 1;
else
	C = CompCladding(Nodes,Trigs,Pdis,'shape');
    dNodesdx = varargin{1};
    dNd = reshape(dNodesdx(:,2:4,:),3*size(Nodes,1),size(dNodesdx,3));
end
fC = C*dNd;
end


function [u] = CompCladding(Nodes,Trigs,Load,varargin)
E =  90e06;   nu = 0.3;  h = 0.02;  prp = [h E nu];
sub = 6; % must be a multiple of 3;
info = SubTriginfoGS(sub);

iD = info{1}; iK = info{2}; iS = info{3}; iL = info{4}; n = info{5};
ID = 6*(5*sub/6*(sub/3+1))+3;

Kall = @(Ktri,elim) sparse([iK.i;iS.i],[iK.j;iS.j],[Ktri(:);iS.v*elim],iD.n,iD.n);
nTri = size(Trigs,1);
if nargin < 4
   u = zeros(nTri,1);
elseif strcmp(varargin{1},'shape')
   u = zeros(nTri,3*size(Nodes,1));
end 
for i = 1:nTri
    NdT  = Trigs(i,2:4);
    Cor  = Nodes(NdT,2:4);
    LdID = Addkron((1:6),6*(NdT-1));
    Pglob = reshape(Load(LdID),3,6);
    AA = Area_Trig(Cor);
    iL.v = iL.f*iL.l*Pglob;
    Lm = AA*iL.v;
    P = sparse(iL.i,iL.j,Lm(:),iD.n,1);
    U = zeros(iD.n,1);
    K = Kall(KtotShell(@(x) Kshell(Cor,prp,n,x),n),0);
    U(iD.dof)  = K(iD.dof,iD.dof)\P(iD.dof);
    ind = 1;
    if nargin < 4    
        out = U(ID);
    elseif strcmp(varargin{1},'shape')
       ind = Addkron(size(Nodes,1)*(0:2),NdT); 
       dKa = @(f,st) dK_all(@(x) Kshell(Cor,prp,n,f,st,x));
       dKp = {dKa(-1,'dKp'),dKa(1,'dKp')};
       dKm = {dKa(-1,'dKm'),dKa(1,'dKm')};
      
       dK = @(d) Kall(KtotShell(@(x) Kshell(Cor,prp,n,x,'shape',d,dKp{1.5+x/2},dKm{1.5+x/2}),n),0);
       dLv = @(d) reshape(Area_Trig(Cor,'shape',d)*iL.v,iD.n,1);
       dP = @(d) sparse(iL.i,iL.j,dLv(d),iD.n,1);
       dV = @(d) dP(d)-dK(d)*U;
       
       daux = [dV(1),dV(2),dV(3),dV(4),dV(5),dV(6),dV(7),dV(8),dV(9)];

       tt  = K(iD.dof,iD.dof)\daux(iD.dof,:);
       out = tt(ID,:);
    end  
  u(i,ind) = u(i,ind)+out;
end
end

function [Ke] = KtotShell(Ke,n)
Rep = @(x,f) repmat(reshape(f(x),324,1),(n+x)*n/2,1);
Ke  = [Rep(1,Ke);Rep(-1,Ke)];
end

function [K] = Kshell(Cor,prp,n,f,varargin)
zm = zeros(1,3);
Tp = @(t) kron(eye(3),[t(3,:),zm;zm,t(2,:);zm,-t(1,:)]);
Tm = @(t) kron(eye(3),[t(1,:),zm;t(2,:),zm;zm,t(3,:)]);
h = prp(1); E = prp(2); nu = prp(3);
locC   = f*Coords_Trig(Cor)/n;
Te = Te_Trig(Cor);

Kp =  KeLCS_TrigPlateDKT(locC(:,1:2),h,E,nu);
Km =  KeLCS_TrigMembraneDrill(locC(:,1:2),h,E,nu);

if nargin < 5
    Kd  = @(Ke,Te)  Te'*Ke*Te;
    K = Kd(Kp,Tp(Te))+Kd(Km,Tm(Te));
elseif strcmp(varargin{1},'dKp')
    K = DiffKeLCS_DKT_Coor(locC,h,E,nu,varargin{2});
elseif strcmp(varargin{1},'dKm')
    K = DiffKeLCS_TMD_Coor(locC,h,E,nu,varargin{2});
elseif strcmp(varargin{1},'shape')
    d = varargin{2};
    dKpa = varargin{3};
    dKma = varargin{4};
    dK = @(Ke,Te,dKe,dTe) dTe'*Ke*Te + Te'*dKe*Te + Te'*Ke*dTe;
    dCs = f*Coords_Trig(Cor,'shape',d)/n;
    dTe = Te_Trig(Cor,'shape',d);
    dKp = dK(Kp,Tp(Te),dKglob(dKpa,dCs),Tp(dTe));
    dKm = dK(Km,Tm(Te),dKglob(dKma,dCs),Tm(dTe));
    K = (dKp+dKm);
end
end

function [dKa] = dK_all(fK)
dKa = cell(9,1);
for i =1:9
    dKa{i}=fK(i);
end
end

function [dK] = dKglob(dKa,dCs)
dK = 0;
for i = 1:9
    dKe = dKa{i};
    dK = dK + dKe*dCs(i);
end
end
