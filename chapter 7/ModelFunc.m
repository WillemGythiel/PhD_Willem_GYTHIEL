function [Model] = ModelFunc(varargin)
if strcmp(varargin{1},'init')
    Model = initiateModel(varargin{2});
elseif strcmp(varargin{1},'update')
    Model = updateModel(varargin{2},varargin{3},varargin{4});
elseif strcmp(varargin{1},'postMerge')
    Model = GeoToModel(varargin{2},varargin{3},varargin{4},varargin{5});
    Model.info = varargin{6};
else
    print('invalid input');
end
end

function [Model] = initiateModel(subDiv)
[Node,Elem,Trig] = GroundStructure(subDiv);
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
Model.Elements  = makeSym(Model,'Bars');
Model.Trigs     = makeSym(Model,'Trigs');
Model.Types     = {1 'beam'};
nSec = size(unique(Model.Elem(:,3)),1);
if size(xSec,1)==1
    xSec = repmat(xSec,nSec,1);
end
if size(xSec,1)~=nSec
    xSec = xSec(end-nSec+1:end,:);
end
Model.Sections  = Model.cache.SecProp((1:nSec)',xSec);
Model.Materials = [1  210e6 0.3  7850 1];
Model.dofs      = makeDOFS(Model.Nodes,Model.Elements,Model.Types);
LC2 = 1;
Model.Pdis      = makeLoad(Model.Nodes,LC2);
Model.cache.get_protected_nodes = @(Node) get_protected_nodes(Node);
end


function [Model] = updateModel(Model,Problem,x)
nVr = size(x,1);
nCvar = size(Problem.SID,1);
nSvar = size(Problem.BID,1);
nNd = size(Model.Nodes,1);
xN = x(1:nCvar);
xS = x(nCvar+1:end);
Model.Pdis      = makeLoad(Model.Nodes,Model.info.LC2);
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
f = 0.15;  b = 0.12;  h = 0.4;
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


function [Pdis] = makeLoad(Nodes,LC2)
if LC2 ~= 0
    lc = 2;
    Pdis = zeros(6*size(Nodes,1),lc);
    indP2 = 6*(Nodes(Nodes(:,2)<=17 & Nodes(:,4)> -100,1)-1)+1;
    Pdis(indP2,2) = LC2;
else
    lc = 1;
    Pdis = zeros(6*size(Nodes,1),lc);
end
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


function [Node,Elem,Trig] = GroundStructure(subDiv)
rounding = 5;
w = 17;
%---------------------------------------------------Definition of the nodes
[Nd,~,~] = PolyMesherAlt(@TrigDomain,subDiv,100);
       num  = size(Nd,1); 
      nID  = (1:num)';                       
       nx  = Nd(:,1);                                  
       ny  = Nd(:,2); 
       nz  = 4*(w^(-4))*(w-nx).*(w+nx).*(w-ny).*(w+ny);
Node      = [nID   nx ny nz;...
             num+1 0 0 -1e3]; 
% -----------------------------------------Definition of the glass elements
Tri   = delaunay(round(Nd(:,1),rounding),round(Nd(:,2),rounding));
Trig = [(1:size(Tri,1))',Tri];    
%------------------------------------------------Definition of the elements
Bar  = unique(sort([Trig(:,2:3);Trig(:,3:4);Trig(:,2:2:4)],2),'rows');
nBr = size(Bar,1); cnt = (1:nBr)'; cte = ones(nBr,1); 
Elem = [cnt cte cnt cte Bar ones(size(Bar,1),1)+num];
end







%% AUXILIAIRY FUNCTIONS %%
function [Node,Element,A] = PolyMesherAlt(Domain,NElem,MaxIter,P)
if ~exist('P','var'), P=PolyMshr_RndPtSet(NElem,Domain); end
NElem = size(P,1);
Tol=5e-6; It=0; Err=1; c=1.5;
BdBox = Domain('BdBox');
Area = 0.5*(BdBox(2)-BdBox(1))*(BdBox(4)-BdBox(3));
Pc = P;
while(It<=MaxIter && Err>Tol)
  Alpha = c*sqrt(Area/NElem);
  P = Pc; %Lloyd's update
  R_P = PolyMshr_Rflct(P,NElem,Domain,Alpha);   %Generate the reflections
  [Node,Element] = voronoin([P;R_P]);           %Construct Voronoi diagram
  [Pc,A] = PolyMshr_CntrdPly(Element,Node,NElem);
  Area = sum(abs(A));
  Err = sqrt(sum((A.^2).*sum((Pc-P).*(Pc-P),2)))*NElem/Area^1.5;
end
[Node,Element] = PolyMshr_ExtrNds(NElem,Node,Element);  %Extract node list
[Node,Element] = PolyMshr_CllpsEdgs(Node,Element,0.1);  %Remove small edges
[Node,Element] = PolyMshr_RsqsNds(Node,Element);        %Reoder Nodes

end
%------------------------------------------------- GENERATE RANDOM POINTSET
function P = PolyMshr_RndPtSet(NElem,Domain)
P=zeros(NElem,2); BdBox=Domain('BdBox'); Ctr=0;
while Ctr<NElem  
  Y(:,1) = (BdBox(2)-BdBox(1))*rand(NElem,1)+BdBox(1);
  Y(:,2) = (BdBox(4)-BdBox(3))*rand(NElem,1)+BdBox(3);
  d = Domain('Dist',Y);
  I = find(d(:,end)<0);                 %Index of seeds inside the domain
  NumAdded = min(NElem-Ctr,length(I));  %Number of seeds that can be added
  P(Ctr+1:Ctr+NumAdded,:) = Y(I(1:NumAdded),:);
  Ctr = Ctr+NumAdded;
end
end
%--------------------------------------------------------- REFLECT POINTSET
function R_P = PolyMshr_Rflct(P,NElem,Domain,Alpha)
eps=1e-8; eta=0.9;
d = Domain('Dist',P);  
NBdrySegs = size(d,2)-1;          %Number of constituent bdry segments
n1 = (Domain('Dist',P+repmat([eps,0],NElem,1))-d)/eps;
n2 = (Domain('Dist',P+repmat([0,eps],NElem,1))-d)/eps;
I = abs(d(:,1:NBdrySegs))<Alpha;  %Logical index of seeds near the bdry
P1 = repmat(P(:,1),1,NBdrySegs);  %[NElem x NBdrySegs] extension of P(:,1)
P2 = repmat(P(:,2),1,NBdrySegs);  %[NElem x NBdrySegs] extension of P(:,2)
R_P(:,1) = P1(I)-2*n1(I).*d(I);  
R_P(:,2) = P2(I)-2*n2(I).*d(I);
d_R_P = Domain('Dist',R_P);
J = abs(d_R_P(:,end))>=eta*abs(d(I)) & d_R_P(:,end)>0;
R_P=R_P(J,:); R_P=unique(R_P,'rows');
end
%---------------------------------------------- COMPUTE CENTROID OF POLYGON
function [Pc,A] = PolyMshr_CntrdPly(Element,Node,NElem)
Pc=zeros(NElem,2); A=zeros(NElem,1);
for el = 1:NElem
  vx=Node(Element{el},1); vy=Node(Element{el},2); nv=length(Element{el}); 
  vxS=vx([2:nv 1]); vyS=vy([2:nv 1]); %Shifted vertices
  temp = vx.*vyS - vy.*vxS;
  A(el) = 0.5*sum(temp);
  Pc(el,:) = 1/(6*A(el,1))*[sum((vx+vxS).*temp),sum((vy+vyS).*temp)];
end
end
%------------------------------------------------------- EXTRACT MESH NODES
function [Node,Element] = PolyMshr_ExtrNds(NElem,Node0,Element0)
map = unique([Element0{1:NElem}]);
cNode = 1:size(Node0,1);
cNode(setdiff(cNode,map)) = max(map);
[Node,Element] = PolyMshr_RbldLists(Node0,Element0(1:NElem),cNode);
end
%----------------------------------------------------- COLLAPSE SMALL EDGES
function [Node0,Element0] = PolyMshr_CllpsEdgs(Node0,Element0,Tol)
while(true)
  cEdge = [];
  for el=1:size(Element0,1)
    if size(Element0{el},2)<4, continue; end  %Cannot collapse triangles
    vx=Node0(Element0{el},1); vy=Node0(Element0{el},2); nv=length(vx);
    beta = atan2(vy-sum(vy)/nv, vx-sum(vx)/nv);
    beta = mod(beta([2:end 1]) -beta,2*pi);
    betaIdeal = 2*pi/size(Element0{el},2);
    Edge = [Element0{el}',Element0{el}([2:end 1])'];
    cEdge = [cEdge; Edge(beta<Tol*betaIdeal,:)];
  end
  if (size(cEdge,1)==0), break; end
  cEdge = unique(sort(cEdge,2),'rows');
  cNode = 1:size(Node0,1);
  for i=1:size(cEdge,1)
    cNode(cEdge(i,2)) = cNode(cEdge(i,1));
  end
  [Node0,Element0] = PolyMshr_RbldLists(Node0,Element0,cNode);
end
end
%--------------------------------------------------------- RESEQUENSE NODES
function [Node,Element] = PolyMshr_RsqsNds(Node0,Element0)
NNode0=size(Node0,1); NElem0=size(Element0,1);
ElemLnght=cellfun(@length,Element0); nn=sum(ElemLnght.^2); 
i=zeros(nn,1); j=zeros(nn,1); s=zeros(nn,1); index=0;
for el = 1:NElem0
  eNode=Element0{el}; ElemSet=index+1:index+ElemLnght(el)^2;
  i(ElemSet) = kron(eNode,ones(ElemLnght(el),1))';
  j(ElemSet) = kron(eNode,ones(1,ElemLnght(el)))';
  s(ElemSet) = 1;
  index = index+ElemLnght(el)^2;
end
K = sparse(i,j,s,NNode0, NNode0);
p = symrcm(K);
cNode(p(1:NNode0))=1:NNode0;
[Node,Element] = PolyMshr_RbldLists(Node0,Element0,cNode);
end
%------------------------------------------------------------ REBUILD LISTS
function [Node,Element] = PolyMshr_RbldLists(Node0,Element0,cNode)
Element = cell(size(Element0,1),1);
[~,ix,jx] = unique(cNode);
if ~isequal(size(jx),size(cNode)), jx=jx'; end % +R2013a compatibility fix
if size(Node0,1)>length(ix), ix(end)=max(cNode); end
Node = Node0(ix,:); 
for el=1:size(Element0,1)
  Element{el} = unique(jx(Element0{el}));
  vx=Node(Element{el},1); vy=Node(Element{el},2); nv=length(vx);
  [~,iix] = sort(atan2(vy-sum(vy)/nv,vx-sum(vx)/nv));
  Element{el} = Element{el}(iix);
end
end
