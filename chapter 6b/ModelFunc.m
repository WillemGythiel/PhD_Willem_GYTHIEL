function [Model] = ModelFunc(varargin)
if strcmp(varargin{1},'init')
    Model = InitiateModel(varargin{2});
elseif strcmp(varargin{1},'update')
    Model = UpdateModel(varargin{2},varargin{3},varargin{4});
else
    print('invalid input');
end
end

function Model = InitiateModel(info)
Model = struct; Model.cache = struct;
pu=4; pv=4; u=(0:info.nSub+1)/(info.nSub+1); v=u;
Model.spline    = @(X) bspline2(reshape(X,info.nCon+1,info.nCon+1),pu,pv,u,v);
Model.NodesCon  = makeControlPoints(info.domain,info.nCon);
Model.Nodes     = ControlPoints2Nodes(Model.NodesCon,info.refNode,info.suppNode,info.loadNode,Model.spline);
Model.Elements  = makeElems(Model.Nodes(:,2:3),info.lMaxF);
Model.Trig      = makeTrigs(Model.Nodes(1:end-1,2:3));
Model.Trigs     = Model.Trig;
Model.cache     = sectionProps(Model.cache);
Model.Types     = {1 'beam'};
Model.Materials = [1  210e6 0.3  7850 1];
Model.dofs      = makeDOFS(Model.Nodes,Model.Elements,Model.Types,info.suppNode,info.suppRule);
Model.Pdis      = makeLoad(Model.Nodes,info.loadNode,info.loadRule);
nSec            = size(unique(Model.Elements(:,3)),1);
Model.Sections  = Model.cache.SecProp((1:nSec)',ones(nSec,1));
Model.cache.get_protected_nodes = @(Node) get_protected_nodes(Node,Model.Elements,info);
Model.info      = info;
Model.info.Rad  = makeRad(Model.Nodes,Model.Elements,info);
end


function Model = UpdateModel(Model,Problem,x)
nVr = size(x,1);
nSvar = size(Problem.BID,1);
xS = x(Problem.nCvar+1:end);
warning('off','all')
dNodesdx = zeros([size(Model.Nodes),nVr]);
if Problem.nCvar>0
    % Update Control points
    Model.NodesCon(Problem.CID)=x(1:Problem.nCvar);
    [Model.Nodes,dNodesdx] = ControlPoints2Nodes(Model.NodesCon,Model.info.refNode,Model.info.suppNode,Model.info.loadNode,Model.spline,Problem);
end
Model.dNodesdx = dNodesdx;
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
warning('on','all')
end



function nodes = makeControlPoints(domain,nCon)
nNodes   = size(domain.Vertices,1);
if nNodes==4
    shapeF = shapeFuncQuad(nCon);
else
    print('Error: only surfaces with quad domain are defined...for now!');
end
nodesXY  = shapeF*domain.Vertices;
nodesZ   = 4*prod(shapeF,2)/max(prod(shapeF,2));
nodesXYZ = [nodesXY,nodesZ];
nodesID  = (1:size(nodesXYZ,1))';
nodes    = [nodesID,nodesXYZ];
end


function [Nodes,varargout] = ControlPoints2Nodes(NodesCon,refNode,suppNode,loadNode,spline,varargin)
[X,dX,~,~] = spline(NodesCon(:,2));
[Y,dY,~,~] = spline(NodesCon(:,3));
[Z,dZ,~,~] = spline(NodesCon(:,4));
nodesXYZ = [X(:),Y(:),Z(:)];
suppNode = cell2mat(suppNode(:,1));
loadNode = cell2mat(loadNode(:,1));
nodesXYZ = unique(round([nodesXYZ;suppNode;loadNode;refNode],5),'rows','stable');
nodesID  = (1:size(nodesXYZ,1))';
Nodes    = [nodesID,nodesXYZ];
if nargin>5 
    Problem = varargin{1};
    dNodesdxAux = zeros([size(Nodes),3*size(dX,2)+size(Problem.BID,1)]);
    [nR,nC] = size(dX);
    dNodesdxAux(1:nR,2,1:3:3*nC) = reshape(dX,nR,1,nC);
    dNodesdxAux(1:nR,3,2:3:3*nC) = reshape(dY,nR,1,nC);
    dNodesdxAux(1:nR,4,3:3:3*nC) = reshape(dZ,nR,1,nC);
    dNodesdx = zeros([size(Nodes),size(Problem.SID,1)+size(Problem.BID,1)]);
    dNodesdx(:,:,1:size(Problem.SID)) = dNodesdxAux(:,:,Problem.SID);
    varargout{1} = dNodesdx;
end
end


function elements = makeElems(nodesXY,lMaxF)
nNodes = size(nodesXY,1);
nodeID = (1:size(nodesXY,1)-1);
elmAll = nchoosek(nodeID,2);
length = sum((nodesXY(elmAll(:,1),:)-nodesXY(elmAll(:,2),:)).^2,2).^0.5;
Lmax   = lMaxF*min(length(length>1e-5));
elemID = elmAll(length<Lmax,:);
nelm   = size(elemID,1);
EltID  = (1:nelm)';
TypID  = ones(nelm,1);
SecID  = ones(nelm,1);
% SecID  = (1:nelm)';
MatID  = ones(nelm,1);
n3     = zeros(nelm,1)+nNodes;
elements = [EltID TypID SecID MatID elemID n3];
end

function trigs = makeTrigs(nodesXY)
trig = delaunay(nodesXY);
trigs = [(1:size(trig))',trig];
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
    pDis(ind,lc) = pDis(ind,lc) + values(:);
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
f = 0.15;  b = 0.12;  h = 0.4;
cache.f = f; cache.b = b; cache.h = h;
APR            = @(x) x*h*b*(1-(1-2*f)^2);
IxPR           = @(x) x*h*b*(1-f)^4/f; 
IyPR           = @(x) x.^3*h^3*b/12*(1-(1-2*f)^4);  
IzPR           = @(x) x*h*b^3/12*(1-(1-2*f)^4);
infM           = @(x) 0*x+Inf;
ZeroX          = @(x) 0*x;
SecProp        = @(ID, x) [ID APR(x) infM(x) infM(x) IxPR(x) IyPR(x) IzPR(x)  x*h/2 x*h/2 ZeroX(x)+b/2 ZeroX(x)+b/2];

dAPR           = @(x) h*b*(1-(1-2*f)^2);
dIxPR          = @(x) h*b*(1-f)^4/f; 
dIyPR          = @(x) x.^2*h^3*b/4*(1-(1-2*f)^4);  
dIzPR          = @(x) h*b^3/12*(1-(1-2*f)^4);
dSecProp       = @(x) [ZeroX(x) dAPR(x) ZeroX(x) ZeroX(x) dIxPR(x) dIyPR(x) dIzPR(x) ZeroX(x)+h/2 ZeroX(x)+h/2 ZeroX(x) ZeroX(x)];
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

