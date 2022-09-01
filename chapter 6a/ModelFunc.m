function [Model] = ModelFunc(varargin)
if strcmp(varargin{1},'init')
    Model = initiateModel(varargin{2});                                    % Model = PretModel('init',subDiv);
elseif strcmp(varargin{1},'update')
    Model = updateModel(varargin{2},varargin{3},varargin{4});              % Model = PretModel('update',Model,Problem,x);
elseif strcmp(varargin{1},'postMerge')
    Model = GeoToModel(varargin{2},varargin{3},varargin{4},varargin{5});   % Model = PretModel('postMerge',Node,Elem,Trig,x);
else
    print('invalid input');
end
end

function [Model] = initiateModel(M)
[Node,Elem,Trig] = GroundStructure(M);
x0S = 1;
Model = GeoToModel(Node,Elem,Trig,x0S);
end


function [Model] = GeoToModel(Node,Elem,Trig,xSec)
Model = struct; Model.cache = struct;
Model.Node      = Node;
Model.Elem      = Elem;
Model.Trig      = Trig;
Model.cache     = sectionProps(Model.cache);
Model.cache     = makeSym(Model,'ABC');
Model.Nodes     = makeSym(Model,'Nodes');
Model.Elements  = makeElems(Model.Nodes,1.5);
Model.Trigs     = makeSym(Model,'Trigs');
Model.Types     = {1 'beam'};
nSec = size(unique(Model.Elem(:,3)),1);
if size(xSec,1)==1
    xSec = repmat(xSec,nSec,1);
end
Model.Sections  = Model.cache.SecProp((1:nSec)',xSec);
Model.Materials = [1  210e6 0.3  7850 1];
Model.dofs      = makeDOFS(Model.Nodes,Model.Elements,Model.Types);
Model.Pdis      = makeLoad(Model.Nodes);
Model.cache.get_protected_nodes = @(Node) get_protected_nodes(Node);
end

function elements = makeElems(Nodes,lMaxF)
nodeID = Nodes(abs(Nodes(:,4))<100,1);
elmAll = nchoosek(nodeID,2);
length = sum((Nodes(elmAll(:,1),2:3)-Nodes(elmAll(:,2),2:3)).^2,2).^0.5;
Lmax   = lMaxF*min(length(length>1e-5));
elemID = elmAll(length<Lmax,:);
nelm   = size(elemID,1);
EltID  = (1:nelm)';
TypID  = ones(nelm,1);
SecID  = ones(nelm,1);
% SecID  = (1:nelm)';
MatID  = ones(nelm,1);
n3     = zeros(nelm,1)+setdiff(Nodes(:,1),nodeID);
elements = [EltID TypID SecID MatID elemID n3];
end


function [Model] = updateModel(Model,Problem,x)
nVr = size(x,1);
nCvar = size(Problem.SID,1);
nSvar = size(Problem.BID,1);
nNd = size(Model.Nodes,1);
xN = x(1:nCvar);
xS = x(nCvar+1:end);
% Update Node(s)
Model.Node(Problem.CID)=xN;
Model.Node(Problem.DID,3) = Model.Node(Problem.DID,2);
Model.Nodes = makeSym(Model,'Nodes');
% generate dNodesdx
dNodesdx = zeros([size(Model.Nodes),nVr]);
checkID  = 3*(Problem.DID-1)+1;
for i = 1:nCvar
    ind = Problem.SID(i);
    dNodesdx(:,2:4,i)=reshape(Model.cache.B(:,ind),3,nNd)';
    if sum(ind==checkID)>0
        dNodesdx(:,2:4,i)=reshape(Model.cache.B(:,ind+1),3,nNd)';
        dNodesdx(:,2:4,i)=dNodesdx(:,2:4,i)+reshape(Model.cache.B(:,ind),3,nNd)';
    end
end
Model.dNodesdx = dNodesdx;
% Update Sections
Model.Sections(Problem.BID,:) = Model.cache.SecProp(Problem.BID,xS);
% Update dSectionsdx
dSectionsdx = zeros([size(Model.Sections),nVr]);
for i=1:nSvar
    dSectionsdx(Problem.BID(i),:,nCvar+i) = Model.cache.dSecProp(xS(i));
end
Model.dSectionsdx = dSectionsdx;
end



function [cache] = sectionProps(cache)
f = 0.15;  b = 0.12;  h = 0.40;
cache.f = f; cache.b = b; cache.h = h;
APR  = @(x) x*h*b*(1-(1-2*f)^2);
IxPR = @(x) x*h*b*(1-f)^4/f; 
IyPR = @(x) x.^3*h^3*b/12*(1-(1-2*f)^4);  
IzPR = @(x) x*h*b^3/12*(1-(1-2*f)^4);
infM = @(x) 0*x+Inf;
ZeroX = @(x) 0*x;
SecProp = @(ID, x) [ID APR(x) infM(x) infM(x) IxPR(x) IyPR(x) IzPR(x)  x*h/2 x*h/2 ZeroX(x)+b/2 ZeroX(x)+b/2];

dAPR  = @(x) h*b*(1-(1-2*f)^2);
dIxPR = @(x) h*b*(1-f)^4/f; 
dIyPR = @(x) x.^2*h^3*b/4*(1-(1-2*f)^4);  
dIzPR = @(x) h*b^3/12*(1-(1-2*f)^4);
dSecProp = @(x) [ZeroX(x) dAPR(x) ZeroX(x) ZeroX(x) dIxPR(x) dIyPR(x) dIzPR(x) ZeroX(x)+h/2 ZeroX(x)+h/2 ZeroX(x) ZeroX(x)];
cache.SecProp  = SecProp;
cache.dSecProp = dSecProp;
end

function [dofs] = makeDOFS(Nodes,Elements,Types)
err = 1e-2;
Nd = Nodes(setdiff(Nodes(:,1),unique(Elements(:,7))),:);
w = (max(Nodes(:,2))-min(Nodes(:,2)))/2;
vw  = @(d) abs(Nd(:,d))>w-err;
NonDofs = [Nd((vw(2)|vw(3)),1)+0.03;
           Nd((vw(2)&vw(3)),1)+0.01;
           Nd((vw(2)&vw(3)),1)+0.02];
dofs = getdof(Elements,Types);
dofs = removedof(dofs,NonDofs);
end


function [Pdis] = makeLoad(Nodes)
lc = 1;
Pdis = zeros(6*size(Nodes,1),lc);
% LC 1
indP1 = 6*(Nodes(Nodes(:,4)> -100,1)-1)+3;   
Pdis(indP1,1)=-2;
end


function [Vertices,Perimeter] = get_protected_nodes(Node)
error = 1e-2;
NdID = Node(1:end-1,1);
Per = @(d) abs(Node(:,d)-17/2)>17/2-error;
interior = NdID(~Per(2)&~Per(3)& abs(Node(:,2)-Node(:,3))>=error);
diameter = NdID(~Per(2)&~Per(3)& abs(Node(:,2)-Node(:,3))< error);
perimX = NdID(Per(2)&~Per(3));
perimY = NdID(Per(3)&~Per(2));
Perimeter = [perimX;perimY;diameter];
Vertices  = setdiff(NdID,[interior;Perimeter]);
end

function [Node,Elem,Trig] = GroundStructure(Model)
rounding = 5;
%---------------------------------------------------Definition of the nodes
Model.Nodes(:,2) = Model.Nodes(:,2)-17;
Model.Nodes(:,3) = Model.Nodes(:,3)-17;
vw1 = @(M,d) abs(M(:,d)-17/2)<=17/2;
test = @(M) M(vw1(M,2)&vw1(M,3)&M(:,2)>=M(:,3),:);
Node = test(Model.Nodes); Node(:,1) = (1:size(Node,1))';
Nd = Node(1:end-1,2:4);
num  = size(Nd,1); 
% -----------------------------------------Definition of the glass elements
Tri   = delaunay(round(Nd(:,1),rounding),round(Nd(:,2),rounding));
Trig = [(1:size(Tri,1))',Tri];    
%------------------------------------------------Definition of the elements
Bar  = unique(sort([Trig(:,2:3);Trig(:,3:4);Trig(:,2:2:4)],2),'rows');
nBr = size(Bar,1); cnt = (1:nBr)'; cte = ones(nBr,1); 
Elem = [cnt cte cnt cte Bar ones(size(Bar,1),1)+num];
end

function [varargout] = makeSym(Model, varargin)
if strcmp(varargin{1},'ABC')
    varargout{1} = makeABC(Model.Node,Model.Elem,Model.cache);      % Model.cache = makeSym(Model,'ABC')
elseif strcmp(varargin{1},'Nodes')
    varargout{1} = mirrorNodes(Model.Node,Model.cache); % Nodes = makeSym(Model,'Nodes')
elseif strcmp(varargin{1},'Trigs')
    varargout{1} = mirrorTrigs(Model.Node,Model.Trig,Model.cache); % Trigs = makeSym(Model,'Trigs')
elseif strcmp(varargin{1},'Bars')
    varargout{1} = mirrorBars(Model.Node,Model.Elem,Model.cache); % Elements = makeSym(Model,'Bars')
else
    print('invalid input');
end
end
 
function [Nodes] = mirrorNodes(Nd,cache)
NodesAux = Nd(:,2:4)';
NodesTT  = cache.B*NodesAux(:);
nNd = size(cache.B,1)/3;
Nodes = zeros(nNd,4);
Nodes(:,1)   = (1:nNd)';
Nodes(:,2:4) = reshape(NodesTT,3,nNd)';
end
 
function [Trigs] = mirrorTrigs(Nd,Trigs,cache)
Tri = Trigs(:,2:4); propTri = Trigs(:,5:end);
TrTT = vertcat(Tri(:,[1,3,2]),Tri); 
Trig = cache.icN(repmat(TrTT,4,1)+kron(size(Nd,1)*(0:7)',ones(size(Tri,1),1)));
Trigs = [(1:size(Trig,1))',Trig,repmat(propTri,8,1)]; 
end
 
function [Elem] = mirrorBars(Nd,Elem,cache)
Bar = Elem(:,5:7); propBar = Elem(:,2:4);
Bar = cache.icN(repmat(Bar,8,1)+kron(size(Nd,1)*(0:7)',ones(size(Bar,1),1)));
BarTT = [repmat(propBar,8,1),Bar];
BarTT = BarTT(cache.iaB,:);
Elem = [(1:size(BarTT,1))',BarTT];
end
 
function [cache] = makeABC(Nd,Elem,cache)
nxa = [Nd(:,2);Nd(:,3)];        nya  = [Nd(:,3);Nd(:,2)];
nxb  = [nxa;-nya;-nxa;nya];     nyb = [nya;nxa;-nya;-nxa];  
rounding = 5;
[~,iaN,icN] = unique(round([nxb, nyb,repmat(Nd(:,4),8,1)],rounding),'rows','stable');
nN = size(Nd,1);
 i1 = (1:8*nN);
  O = ones(2*size(Nd,1),1);
 iC = 3*(0:size(Nd,1)-1);
 iX = iC+1; iY = iC+2; iZ = iC+3; 
 jx = [iX,iY,iY,iX,iX,iY,iY,iX];
 jy = [iY,iX,iX,iY,iY,iX,iX,iY];
 jz = [iZ,iZ,iZ,iZ,iZ,iZ,iZ,iZ];
 v1x = [O;-O;-O;O];
 v1y = [O;O;-O;-O];
 v1z = [O;O;O;O];
makeB = @(jI,vI) sparse(i1(:),jI(:),vI(:),size(i1(:),1),3*nN);
Bx = makeB(jx,v1x); Bx = Bx(iaN,:);
By = makeB(jy,v1y); By = By(iaN,:);
Bz = makeB(jz,v1z); Bz = Bz(iaN,:);
B = sparse(3*size(Bx,1),size(Bx,2)); 
B(1:3:end,:)=Bx; B(2:3:end,:)=By; B(3:3:end,:)= Bz;
cache.B = B; 
cache.indsN = icN(1:size(Nd,1)); cache.icN = icN;
Bar = Elem(:,5:7); 
Bar = cache.icN(repmat(Bar,8,1)+kron(size(Nd,1)*(0:7)',ones(size(Bar,1),1)));
[~,iaB,~] = unique(Bar,'rows','stable');
ci = (1:size(iaB,1));
cjTT = repmat((1:size(Elem,1))',8,1); cj = cjTT(iaB);
cv = ones(size(iaB,1),1);
cache.C = sparse(ci,cj,cv,size(iaB,1),size(Elem,1));
cache.iaB = iaB;
end
