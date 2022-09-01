function [Model] = ModelFunc(varargin)
if strcmp(varargin{1},'init')
    Model = InitiateModel(varargin{2});                                    
elseif strcmp(varargin{1},'update')
    Model = UpdateModel(varargin{2},varargin{3},varargin{4});             
elseif strcmp(varargin{1},'topology')
    Model = TopolgyModel(varargin{2},varargin{3},varargin{4});  
elseif strcmp(varargin{1},'changeLoad')
    Model = ChangeLoad(varargin{2},varargin{3});  
else
    print('invalid input');
end
end

function Model = InitiateModel(info)
Model = struct; Model.cache = struct;
Model.Node      = makeNodes(info.domain,info.refNode,info.suppNode,info.loadNode,info.nStep);
[Model.Nodes,Model.link] = makeSym(Model.Node);
Model.Elements  = makeElems(Model.Nodes(:,2:3),info.lMaxF,Model.link);
Model.Trig      = makeTrigs(Model.Nodes(1:end-1,2:3));
Model.Trigs     = Model.Trig;
Model.cache     = sectionProps(Model.cache);
Model.Types     = {1 'beam'};
Model.Materials = [1  210e6 0.3  7850 1];
Model.dofs      = makeDOFS(Model.Nodes,Model.Elements,Model.Types,info.suppNode,info.suppRule);
Model.Pdis      = makeLoad(Model.Nodes,info.loadNode,info.loadRule);
nSec            = size(unique(Model.Elements(:,3)),1);
Model.Sections  = Model.cache.SecProp((1:nSec)',ones(nSec,1)*.5);
Model.cache.get_protected_nodes = @(Node) get_protected_nodes(Node,Model.Elements,info);
Model.info      = info;
Model.info.Rad  = makeRad(Model.Node(1:end-1,:),Model.Elements,info);
end

function Model = ChangeLoad(Model,info)
Model.Pdis      = makeLoad(Model.Nodes,info.loadNode,info.loadRule);
Model.Types     = {1 'beam'};
Model.Materials = [1  210e6 0.3  7850 1];
nSec            = size(unique(Model.Elements(:,3)),1);
Model.cache     = sectionProps(Model.cache);
Model.Sections  = Model.cache.SecProp((1:nSec)',ones(nSec,1)*.5);
end



function Model = UpdateModel(Model,Problem,x)
nVr = size(x,1);
nSvar = size(Problem.BID,1);
xS = x(Problem.nCvar+1:end);
% generate dNodesdx
Model.dNode = Problem.dNode;
if Problem.nCvar>0
    % Update Node(s)
    Aux = zeros(size(Model.Node));
    for i=1:nVr
        Aux = Aux+Model.dNode(:,:,i)*x(i);
    end
    Model.Node=Problem.NodL+Aux;
end
[Model.Nodes,Model.link,Model.dNodesdx] = makeSym(Model.Node,Model.dNode);
% Update dSectionsdx
dSectionsdx = zeros([size(Model.Sections),nVr]);
if nSvar>0
    % Update Sections
    Model.Sections(Problem.BID,:) = Model.cache.SecProp(Problem.BID,xS);
    for i=1:nSvar
        dSectionsdx(Problem.BID(i),:,Problem.nCvar+i) = Model.cache.dSecProp(xS(i));
    end
end
Model.dSectionsdx = dSectionsdx;
Model.Trigs     = Model.Trig;
end

function [Model] = TopolgyModel(Model,Problem,x)
[Model.Sections,Model.Elements] = Topology_Selection(Model,Problem,x);
end


function nodes = makeNodes(domain,refNode,suppNode,loadNode,nStep)
nNodes   = size(domain.Vertices,1);
if nNodes==4
    shapeF = shapeFuncQuad(nStep);
else
    subDiv   = (0:nStep)/nStep;
    shapeF   = unique(nchoosek(repmat(subDiv,1,nNodes),nNodes),'rows');
    shapeF   = shapeF(abs(sum(shapeF,2)-1)<1e-5,:);
end

shapeF = shapeF+repmat(prod(shapeF,2)~=0,1,3).*(rand(size(shapeF))-0.5)*0.20;
nodesXY  = shapeF*domain.Vertices;
suppNode = cell2mat(suppNode(:,1));
loadNode = cell2mat(loadNode(:,1));
nodesZ   = 0*nodesXY(:,1);
nodesXYZ = [nodesXY,nodesZ];
nodesXYZ = unique(round([nodesXYZ;suppNode;loadNode;refNode],5),'rows','stable');
nodesID  = (1:size(nodesXYZ,1))';
nodes    = [nodesID,nodesXYZ];

[f,~] = surfaceSphere(nodes,domain.Vertices);
nodes(1:end-1,4) = f(1:end-1);
end

function elements = makeElems(nodesXY,lMaxF,link)
nNodes = size(nodesXY,1);
nodeID = (1:size(nodesXY,1)-1);
Cn = nchoosek(nodeID,2);
n1  = Cn(:,1);
n2  = Cn(:,2);
dc  = nodesXY(n2,:)-nodesXY(n1,:);
L  = sum(dc.^2,2).^0.5;
Cs = dc./repmat(L,1,2);
Lmax = lMaxF;
keep = L<Lmax;
for i=1:size(Cn)-1
    if keep(i)                                                             % only check if you intend to keep the first bar
        for j=i+1:size(Cn)
            if numel(nonzeros(Cn(i,:)==Cn(j,:)))>0                         % only check if both bars share a node
                if keep(j)                                                 % only check if you intend to keep the second bar
                    if sum(abs(Cs(i,:)-Cs(j,:)),2)<1e-2                         % verify if parallel
                        if keep(i) %doublecheck
                            keep(i) = L(i)<L(j);                               % keep bar i if it is shorter
                        end
                        keep(j) = L(j)<L(i);                               % keep bar j if it is shorter
                    end
                end
            end
        end
    end
end
elemID = Cn(keep==true,:);

IDaux = zeros(size(elemID));
for k=1:size(link,1)
    lk = link{k};
    if size(lk,1)>0
        for l=1:size(lk,1)
            IDaux(abs(elemID-lk(l))<1e-2)=k;
        end
    end
end
[~,~,SecID]=unique([sort(IDaux,2),round(L(keep==true),2)],'rows');
nelm   = size(elemID,1);
EltID  = (1:nelm)';
TypID  = ones(nelm,1);
MatID  = ones(nelm,1);
n3     = zeros(nelm,1)+nNodes;
elements = [EltID TypID SecID MatID elemID n3];
end

function trigs = makeTrigs(nodesXY)
trig = delaunay(round(nodesXY,5));
cs = @(dim,nd) nodesXY(trig(:,nd),dim);
A = ((cs(1,1)-cs(1,2)).*(cs(2,1)-cs(2,3))-(cs(1,1)-cs(1,3)).*(cs(2,1)-cs(2,2)));
trig = trig(A.^2>1e-10,:);
trigs = [(1:size(trig,1))',trig];
end

function dofs = makeDOFS(Nodes,Elements,Types,suppNode,suppRule)
nSupp   = size(suppNode,1); nRule = size(suppRule,1); 
suppDOF = cell(nSupp+nRule,1);
eq      = @(n1,n2) abs(n1(:,1)-n2(:,1))+abs(n1(:,2)-n2(:,2))<1e-5;
for i=1:nSupp
    suppID     = Nodes(eq(Nodes(:,2:3),suppNode{i,1}),1);
    suppDOF{i} = suppID+suppNode{i,2}';
end
for i=1:nRule
    rule             = suppRule{i,1};
    suppID           = rule(Nodes);
    suppDOF{i+nSupp} = Addkron(suppID,suppRule{i,2}');
end
NonDofs = cell2mat(suppDOF);
allDofs = getdof(Elements,Types);
dofs    = removedof(allDofs,NonDofs);
end


function pDis = makeLoad(Nodes,loadNode,loadRule)
nPts  = size(loadNode,1); nRule = size(loadRule,1); 
nLC   = max([cell2mat(loadNode(:,3));cell2mat(loadRule(:,3))]);
pDis  = zeros(6*size(Nodes,1),nLC);
eq    = @(n1,n2) abs(n1(:,1)-n2(:,1))+abs(n1(:,2)-n2(:,2))<1e-5;
Index = @(id) Addkron(6*(id-1),(1:6)');
for i=1:nPts
    loadID       = Nodes(eq(Nodes(:,2:3),loadNode{i,1}),1);
    values       = loadNode{i,2};
    lc           = loadNode{i,3};
    ind          = Index(loadID);
    if numel(ind)>0
        pDis(ind,lc) = pDis(ind,lc) + values(:);
    end
end
for i=1:nRule
    Rule         = loadRule{i,1};
    values       = loadRule{i,2};
    lc           = loadRule{i,3};
    loadID       = Rule(Nodes);
    ind          = Index(loadID);
    pDis(ind,lc) = pDis(ind,lc) + repmat(values(:),size(loadID,1),1);
end
end

function [cache] = sectionProps(cache)
f = 0.10;  Rmax = 0.3;
cache.f = f; cache.b = Rmax/2; cache.h = Rmax/2; cache.R = Rmax;
APR            = @(x) pi*x.^2*Rmax^2*(1-(1-f)^2);
IxPR           = @(x) pi*x.^4*Rmax^4/2*(1-(1-f)^4);
IyPR           = @(x) pi*x.^4*Rmax^4/4*(1-(1-f)^4);  
IzPR           = @(x) pi*x.^4*Rmax^4/4*(1-(1-f)^4);
infM           = @(x) 0*x+Inf;
ZeroX          = @(x) 0*x;
SecProp        = @(ID, x) [ID APR(x) infM(x) infM(x) IxPR(x) IyPR(x) IzPR(x)  x*Rmax x*Rmax ZeroX(x)+Rmax ZeroX(x)+Rmax];

dAPR           = @(x) pi*2*x*Rmax^2*(1-(1-f)^2);
dIxPR          = @(x) pi*2*x.^3*Rmax^4*(1-(1-f)^4);
dIyPR          = @(x) pi*x.^3*Rmax^4*(1-(1-f)^4);
dIzPR          = @(x) pi*x.^3*Rmax^4*(1-(1-f)^4);
dSecProp       = @(x) [ZeroX(x) dAPR(x) ZeroX(x) ZeroX(x) dIxPR(x) dIyPR(x) dIzPR(x) ZeroX(x)+Rmax ZeroX(x)+Rmax ZeroX(x) ZeroX(x)];
cache.SecProp  = SecProp;
cache.dSecProp = dSecProp;
end

function [rad] = makeRad(Nodes,Elements,info)
Perimeter = get_perimeter(Nodes,Elements,info.domain);
Verts     = unique(cell2mat(Perimeter(:,1)));
circs     = nchoosek(Verts,3);
rad       = max(CircumRadiusTrigGS(Nodes,[circs(:,1),circs]));
end


function [Vertices,Perimeter] = get_protected_nodes(Nodes,Elements,info)
nSupp     = size(info.suppNode,1);
Vertices  = zeros(nSupp,1);
eq        = @(n1,n2) abs(n1(:,1)-n2(:,1))+abs(n1(:,2)-n2(:,2))<1e-5;
for i=1:nSupp
    Vertices(i) = Nodes(eq(Nodes(:,2:3),info.suppNode{i,1}),1);
end
Perimeter = get_perimeter(Nodes,Elements,info.domain);
Perimeter = nonzeros(cell2mat(Perimeter(:,2)));
end

function [Sections,Elements] = Topology_Selection(Model,Problem,x)
xS = x(Problem.nCvar+1:end);
lb = Problem.lb(Problem.nCvar+1:end);
eq = @(x,y) abs(x-y)<1e-4;
eq2 = @(x,y) (x-y)./y<0.2;
sec = Model.Sections(~eq2(xS,lb),:);
elm = cell(size(sec,1),1);
Sections = sec;
for i=1:size(sec,1)
    Sections(i,1) = i;
    elI = Model.Elements(eq(Model.Elements(:,3),sec(i,1)),:);
    elI(:,3) = i;
    elm{i}=elI;
end
Elements = cell2mat(elm);
end

function [Nodes,link,varargout] = makeSym(node,varargin)
x = node(:,2); y = node(:,3); z = node(:,4);
RT =  [(x.^2+y.^2).^0.5,atan2(y,x)];
n0 = size(RT,1);
T1 = min(RT(:,2));
T2 = max(RT(:,2));
delta = round(pi/(T2-T1),0);

Nodes = zeros(2*delta*n0+1,4);
dNd = zeros(2*delta*n0+1,3*n0,3);
sel = zeros(2*delta*n0+1,1);
for i=1:n0-1
    Ri = RT(i,1);
    Ti = RT(i,2);
    if Ri~=min(RT(:,1))
        for d = 0:delta-1
            [Crot,dCrot] = make_Crot(x(i),y(i),z(i),T1,T2,d);
            Nodes(i+d*n0,2:4) = Crot;
            dNd(i+d*n0,(0:2)*n0+i,1:3) = reshape(dCrot',1,3,3);
            sel(i+d*n0) = i+d*n0;
            if abs(Ti-T1)>1e-5&&abs(Ti-T2)>1e-5
                [Cmir,dCmir] = make_Cmir(x(i),y(i),z(i),T1,T2,d);
                Nodes(i+delta*n0+d*n0,2:4) = Cmir;
                sel(i+delta*n0+d*n0) = i+delta*n0+d*n0;
                dNd(i+delta*n0+d*n0,(0:2)*n0+i,1:3) = reshape(dCmir(:),1,3,3);
            end
        end
    else
        if z(i)>-1e-3
            Nodes(i,2:4) = node(i,2:4);
            sel(i)=i;
        end
    end
end
Nodes(end,2:4) = node(n0,2:4);
sel(end)=2*delta*n0+1;
sel = nonzeros(sel);
Nodes = Nodes(sel,:); Nodes(:,1)=(1:size(Nodes,1));
dNd = dNd(sel,:,:);

aux = (1:size(dNd,1));
link = cell(n0,1);
for k = 1:n0
    link{k}=nonzeros(aux(abs(sum(dNd(:,k,:),3))>1e-4));
end

if nargin>1
    dnode = varargin{1};
    dAux = reshape(dnode(:,2:4,:),3*size(dnode,1),size(dnode,3));
    dNodesdx = zeros(size(Nodes,1),4,size(dnode,3));
    for i = 1:3
        dNodesdx(:,i+1,:)=reshape(dNd(:,:,i)*dAux,size(Nodes,1),1,size(dnode,3));
    end
    varargout{1} = dNodesdx;
end
end

function [Cmir,dCmir] = make_Cmir(xi,yi,zi,T1,T2,n)
 Cmir = [ cos(2*n*(T1 - T2) + atan2(yi, xi) - 2*T2)*(xi^2 + yi^2)^(1/2);
         -sin(2*n*(T1 - T2) + atan2(yi, xi) - 2*T2)*(xi^2 + yi^2)^(1/2);
          zi];
dCmir = [  (xi*cos(2*n*(T1 - T2) + atan2(yi, xi) - 2*T2))/(xi^2 + yi^2)^(1/2) + (sin(2*n*(T1 - T2) + atan2(yi, xi) - 2*T2)*(xi^2 + yi^2)^(1/2)*(imag(xi) + real(yi)))/((imag(xi) + real(yi))^2 + (imag(yi) - real(xi))^2),   (yi*cos(2*n*(T1 - T2) + atan2(yi, xi) - 2*T2))/(xi^2 + yi^2)^(1/2) + (sin(2*n*(T1 - T2) + atan2(yi, xi) - 2*T2)*(xi^2 + yi^2)^(1/2)*(imag(yi) - real(xi)))/((imag(xi) + real(yi))^2 + (imag(yi) - real(xi))^2), 0;
         - (xi*sin(2*n*(T1 - T2) + atan2(yi, xi) - 2*T2))/(xi^2 + yi^2)^(1/2) + (cos(2*n*(T1 - T2) + atan2(yi, xi) - 2*T2)*(xi^2 + yi^2)^(1/2)*(imag(xi) + real(yi)))/((imag(xi) + real(yi))^2 + (imag(yi) - real(xi))^2), - (yi*sin(2*n*(T1 - T2) + atan2(yi, xi) - 2*T2))/(xi^2 + yi^2)^(1/2) + (cos(2*n*(T1 - T2) + atan2(yi, xi) - 2*T2)*(xi^2 + yi^2)^(1/2)*(imag(yi) - real(xi)))/((imag(xi) + real(yi))^2 + (imag(yi) - real(xi))^2), 0;
                                                                                                                                                                                                                        0,                                                                                                                                                                                                                0, 1];
end

function [Crot,dCrot] = make_Crot(xi,yi,zi,T1,T2,n)
 Crot = [  cos(atan2(yi, xi) - 2*n*(T1 - T2))*(xi^2 + yi^2)^(1/2);
           sin(atan2(yi, xi) - 2*n*(T1 - T2))*(xi^2 + yi^2)^(1/2);
                                                    zi];
dCrot = [(xi*cos(atan2(yi, xi) - 2*n*(T1 - T2)))/(xi^2 + yi^2)^(1/2) + (sin(atan2(yi, xi) - 2*n*(T1 - T2))*(xi^2 + yi^2)^(1/2)*(imag(xi) + real(yi)))/((imag(xi) + real(yi))^2 + (imag(yi) - real(xi))^2), (yi*cos(atan2(yi, xi) - 2*n*(T1 - T2)))/(xi^2 + yi^2)^(1/2) + (sin(atan2(yi, xi) - 2*n*(T1 - T2))*(xi^2 + yi^2)^(1/2)*(imag(yi) - real(xi)))/((imag(xi) + real(yi))^2 + (imag(yi) - real(xi))^2), 0;
         (xi*sin(atan2(yi, xi) - 2*n*(T1 - T2)))/(xi^2 + yi^2)^(1/2) - (cos(atan2(yi, xi) - 2*n*(T1 - T2))*(xi^2 + yi^2)^(1/2)*(imag(xi) + real(yi)))/((imag(xi) + real(yi))^2 + (imag(yi) - real(xi))^2), (yi*sin(atan2(yi, xi) - 2*n*(T1 - T2)))/(xi^2 + yi^2)^(1/2) - (cos(atan2(yi, xi) - 2*n*(T1 - T2))*(xi^2 + yi^2)^(1/2)*(imag(yi) - real(xi)))/((imag(xi) + real(yi))^2 + (imag(yi) - real(xi))^2), 0;
                                                                                                                                                                                                        0,                                                                                                                                                                                                0, 1];
end

