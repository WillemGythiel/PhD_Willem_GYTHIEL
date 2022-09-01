function [P] = LoadVoronoi(Nodes,Elements,Load,varargin)
Nd = Nodes(setdiff(Nodes(:,1),unique(Elements(:,7))),:);
Tri = delaunay(Nd(:,2),Nd(:,3));
Trigs = [(1:size(Tri))',Nd(Tri(:,1),1),Nd(Tri(:,2),1),Nd(Tri(:,3),1)];%Nd(Tri(:,1),1) to get right node nr becasue some node nrs were removed by previous command
if nargin<4
    Lfunc = @(Cor,in) loadFunc(Cor,in,varargin);
    Ld_Nd = Load_Builder(Nodes,Trigs,Lfunc);
else
    dNodesdx = varargin{1};
    Lfunc = @(Cor,in) loadFunc(Cor,in,'derivative');
    dNd = reshape(dNodesdx(:,2:4,:),3*size(Nodes,1),size(dNodesdx,3));
    Ld_Nd = Load_Builder(Nodes,Trigs,Lfunc,'derivative')*dNd;
end

nLc = size(Load,2);
P = zeros([size(Ld_Nd),nLc]);
for lc = 1:nLc
    Ld_Mag = diag(Load(:,lc));
    P(:,:,lc) = Ld_Mag*Ld_Nd;
end
end

function [P] = Load_Builder(Nodes,Trigs,Lfunc,varargin)
in = center_inside(Nodes,Trigs);
nNd = size(Nodes,1);    nLds = 6*nNd;
if nargin < 4
    nC = 1;
elseif strcmp(varargin{1},'derivative')
    nC = 3*nNd;
else
    print('Input not supported')
end
Pv = zeros(nLds,nC);
for i = 1:size(Trigs,1)
    NdT  =  Trigs(i,2:4);
    Cor  = Nodes(NdT,2:4);
    Li   = Addkron(6*(NdT-1),(1:6));
    if nargin < 4
        Lj = 1;
    else
        Lj = Addkron(NdT,nNd*(0:1));
    end
    Lpn = kron(Lfunc(Cor,in(i)),ones(6,1));
    Pv(Li,Lj) = Pv(Li,Lj)+Lpn;
end
[Pi,Pj] = ndgrid((1:nLds),(1:nC));
P = sparse(Pi(:),Pj(:),Pv(:),nLds,nC);
end

function [A] = loadFunc(Cor,in,varargin)
XYmat = reshape(Cor(:,1:2)',6,1);
A = AreaTriVoronoi(XYmat,in);
if strcmp(varargin{1},'derivative')
    A = AreaTriVoronoi(XYmat,in,'derivative');
end
end

function [in] = center_inside(Nodes,Trigs)
C = @(t,c) Nodes(Trigs(:,t+1),c+1);
x1=C(1,1); y1=C(1,2); x2=C(2,1); y2=C(2,2); x3=C(3,1); y3=C(3,2);

xc=(x1.^2.*y2-x1.^2.*y3-x2.^2.*y1+x2.^2.*y3+x3.^2.*y1-x3.^2.*y2+y1.^2.*y2-y1.^2.*y3-y1.*y2.^2+y1.*y3.^2+y2.^2.*y3-y2.*y3.^2)./(2.*(x1.*y2-x2.*y1-x1.*y3+x3.*y1+x2.*y3-x3.*y2));
yc=(-x1.^2.*x2+x1.^2.*x3+x1.*x2.^2-x1.*x3.^2+x1.*y2.^2-x1.*y3.^2-x2.^2.*x3+x2.*x3.^2-x2.*y1.^2+x2.*y3.^2+x3.*y1.^2-x3.*y2.^2)./(2.*(x1.*y2-x2.*y1-x1.*y3+x3.*y1+x2.*y3-x3.*y2));

k = convhull(Nodes(:,2),Nodes(:,3));
in = inpolygon(xc,yc,Nodes(k,2),Nodes(k,3));
end