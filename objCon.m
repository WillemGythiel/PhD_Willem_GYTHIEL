function [varargout] = objCon(x,ObjCon,func,Mfunc,Lfunc,M,Problem,derivatives,varargin)
x = x(:);
M = Mfunc('update',M,Problem,x);





expectedObjCon = {'objective','constraint'};

p = inputParser; 
addRequired(p,'func');
addRequired(p,'ObjCon',@(x) any(validatestring(x,expectedObjCon)));
parse(p,func,ObjCon);
 
nlc = size(M.Pdis,2);       nVar = size(x,1);           nDof = size(M.dofs,1); 
f0s = cell(14,1);           df0s = cell(14,1);          nElm = size(M.Elements,1);
f0CELL = cell(nlc,1);       df0CELL = cell(nlc,1);
 
dFrac = @(t,dt,n,dn,nVar) (dt.*repmat(n,1,nVar)-repmat(t,1,nVar).*dn)./repmat(n.^2,1,nVar);
 prop = @(col) M.Sections(M.Elements(:,3),col);
dprop = @(col) reshape(M.dSectionsdx(M.Elements(:,3),col,:),nElm,nVar);

if ~isempty(intersect(p.Results.func,{'stress','buckling','complianceCladdingPA',...
                                      'displacements','compliance','Eulerbuckling','loadLim'}))
    [K,~,dKdx] = asmkm(M.Nodes,M.Elements,M.Types,M.Sections,M.Materials,M.dofs,M.dNodesdx,M.dSectionsdx);
    Kinv = K\speye(size(K));
    P = Lfunc(M.Nodes,M.Elements,M.Trigs,M.Pdis);
    DOFs = round(6*(round(M.dofs)-1)+100*mod(M.dofs,1));
    P = reshape(P(DOFs,:,:),nDof,nlc);
    U = Kinv*P;
    if derivatives
        dP = Lfunc(M.Nodes,M.Elements,M.Trigs,M.Pdis,M.dNodesdx);
        dUdx = zeros(size(U,1),nlc,nVar);
        for i = 1:size(x,1)
            dFdx = -dKdx{i}*U+reshape(dP(DOFs,i,:),nDof,nlc);
            dUdx(:,:,i) = Kinv*dFdx;
        end
    end

end

if any(strcmp(p.Results.func,'volume'))
    try 
        Vmax = Problem.Vmax;
    catch
        Vmax = 1;
    end
    pw = M.info.power;
    nVr = size(x,1);
    nSvar = size(Problem.BID,1);
    xS = x(Problem.nCvar+1:end);
    SectionsP(Problem.BID,:) = M.cache.SecProp(Problem.BID,xS.^pw);
    dSectionsP = zeros([size(M.Sections),nVr]);
    if nSvar>0
        for i=1:nSvar
            dSectionsP(Problem.BID(i),:,Problem.nCvar+i) = M.cache.dSecProp(xS(i).^pw)*pw*xS(i)^(pw-1);
        end
    end
    [Ve,dVedx] = elemvolumes(M.Nodes,M.Elements,M.Types,SectionsP,M.dNodesdx,dSectionsP);
    A = prop(2); dA    = dprop(2);
    jointCost = M.info.jointCost;
    f0s{1} = sum(Ve+jointCost*A)/Vmax;
    if derivatives
        df0s{1} = sum(dVedx+jointCost*dA)/Vmax;
    end
end
 
 
if any(strcmp(p.Results.func,'area'))
    f0s{2} = -1e-3*AreaFunc(M.Nodes,M.Trig)+1;
    if derivatives
        df0s{2} = -1e-3*AreaFunc(M.Nodes,M.Trig,M.dNodesdx);
    end
end
 
if ~isempty(intersect(p.Results.func,{'stress','buckling','Eulerbuckling'}))
    if derivatives
        [ForcesLCS,~,dForcesLCS,~] = elemforces(M.Nodes,M.Elements,M.Types,M.Sections,M.Materials,M.dofs,U,[],M.dNodesdx,M.dSectionsdx,dUdx);
    else
        ForcesLCS = elemforces(M.Nodes,M.Elements,M.Types,M.Sections,M.Materials,M.dofs,U,[]);
    end
    
end
 
for lc = 1:nlc
    if ~isempty(intersect(p.Results.func,{'stress','buckling','Eulerbuckling'}))
        frc  = @(col) ForcesLCS(:,col,lc);   
        N = frc(1); My1 = frc(5); My2 = frc(11); Mz1 = frc(6); Mz2 = frc(12);
        if derivatives
            dfrc  = @(col) reshape(dForcesLCS(:,col,lc,:),size(N,1),nVar);
            dN = dfrc(1); dMy1 = dfrc(5); dMy2 = dfrc(11); dMz1 = dfrc(6); dMz2 = dfrc(12);
        end
    end
    %_________________________________________________________displacements
    if any(strcmp(p.Results.func,'displacements'))
        w = M.info.Rad;
        [~,IDu] = selectdof(M.dofs,Problem.int+0.03);
        IDu = nonzeros(IDu);
        nIDu = size(IDu,1);
        fdisp = @(Ux) [Ux;-Ux]*(200/w);
        f0s{3} = fdisp(U(IDu,lc));
        if derivatives
            df0s{3} =fdisp(reshape(dUdx(IDu,lc,:),nIDu,nVar));
        end
    end 
    %______________________________________________________________buckling
    if any(strcmp(p.Results.func,'buckling'))
        [L,dL] = elemsizes(M.Nodes,M.Elements,M.Types,M.dNodesdx);
        Matprop = @(col) M.Materials(M.Elements(:,4),col);
        E = Matprop(2); fy = Problem.SigmaT;
        A = prop(2); Iy = prop(6); Iz = prop(7);  y = prop(10); z = prop(8);
        dA    = dprop(2); dIy = dprop(6); dIz = dprop(7); dy = dprop(10); dz = dprop(8);
        
        [UR,dUR] = bucklingEC(N,My1,My2,Mz1,Mz2,E,A,Iy,Iz,L,y,z,fy,dN,dMy1,dMy2,dMz1,dMz2,dA,dIy,dIz,dL,dy,dz);
        func = UR;
        [fCluster,dfCluster] = clusterFunc(func,M.info.nClusterR,M.info.shift);
        f0s{4} = fCluster;
        if derivatives

            df0s{4} = dfCluster*real(dUR);
        end
    end    
    %________________________________________________________________stress
    if any(strcmp(p.Results.func,'stress'))
                        A = prop(2); Iy = prop(6); Iz = prop(7);  y = prop(10); z = prop(8);
        limitT = Problem.SigmaT;
        NRd  = A*limitT;
        MyRd = limitT*Iy./z;
        MzRd = limitT*Iz./y;
        [I,J,K] = ndgrid([-1,1],[-1,1],[-1,1]);
        UrFunc = @(UrN,UrMy,UrMz) kron(I(:),UrN)+kron(J(:),UrMy)+kron(K(:),UrMz);
        func = [UrFunc(N./NRd,My1./MyRd,Mz1./MzRd);UrFunc(N./NRd,My2./MyRd,Mz2./MzRd)];
       [fCluster,dfCluster] = clusterFunc(func,M.info.nClusterS,M.info.shift);
        f0s{5} = fCluster;
        if derivatives
            dA    = dprop(2); dIy = dprop(6); dIz = dprop(7); dy = dprop(10); dz = dprop(8);
            dNRd  = dA*limitT;
            dMyRd = limitT*dFrac(Iy,dIy,z,dz,nVar);
            dMzRd = limitT*dFrac(Iz,dIz,y,dy,nVar);
            dUrN   = dFrac(N,dN,NRd,dNRd,nVar);  
            dUrMy1 = dFrac(My1,dMy1,MyRd,dMyRd,nVar);
            dUrMy2 = dFrac(My2,dMy2,MyRd,dMyRd,nVar);
            dUrMz1 = dFrac(Mz1,dMz1,MzRd,dMzRd,nVar);
            dUrMz2 = dFrac(Mz2,dMz2,MzRd,dMzRd,nVar);
            dfunc  = [UrFunc(dUrN,dUrMy1,dUrMz1);UrFunc(dUrN,dUrMy2,dUrMz2)];
            
            df0s{5} = dfCluster*dfunc;
        end
    end     
    %____________________________________________________________compliance
    if any(strcmp(p.Results.func,'compliance'))
        try 
            Cmax = Problem.Cmax;
        catch
            Cmax = 1;
        end
        f0s{6} = P(:,lc)'*U(:,lc)/Cmax;
        if derivatives
            dCdx = zeros(1,nVar);
            for k = 1:nVar
                dCdx(k) = -U(:,lc)'*(dKdx{k}*U(:,lc)-dP(DOFs,k,lc));
            end
            df0s{6} = dCdx/Cmax;
        end  
    end  
    %____________________________________________________complianceCladding
    if any(strcmp(p.Results.func,'complianceCladding'))
        try 
            cCmaxC = Problem.cCmaxC;
        catch
            cCmaxC = 1;
        end
        f0s{7} = Compliance_CladdingGS(M.Nodes,M.Trig,M.Pdis)/cCmaxC;
        if derivatives
            df0s{7} = Compliance_CladdingGS(M.Nodes,M.Trig,M.Pdis,M.dNodesdx)/cCmaxC;
        end      
    end 
    %____________________________________________________deflectionCladding
    if any(strcmp(p.Results.func,'deflectionCladding'))
        factor = Problem.factor;
        t = Deflection_CladdingGS(M.Nodes,M.Trig,M.Pdis);
        n = CircumRadiusTrigGS(M.Nodes,M.Trig);
        f0s{8} = factor*t./n;
        if derivatives
            dt = Deflection_CladdingGS(M.Nodes,M.Trig,M.Pdis,M.dNodesdx);
            dn = CircumRadiusTrigGS(M.Nodes,M.Trig,M.dNodesdx);
            df0s{8} = factor*(dt.*repmat(n,1,size(dt,2))-repmat(t,1,size(dn,2)).*dn)./repmat(n,1,size(dt,2)).^2;
        end      
    end 
    %_______________________________________________________surfaceFunction
    if any(strcmp(p.Results.func,'surfaceFunc'))
        errMax = M.info.delta^2;
        nID = M.Node(1:end-1,1);
        f0s{9} = SurfaceFunc(M.Nodes,nID)/errMax;
        if derivatives
            df0s{9} = SurfaceFunc(M.Nodes,nID,M.dNodesdx)/errMax;
        end      
    end
    %__________________________________________________complianceCladdingPA
    if any(strcmp(p.Results.func,'complianceCladdingPA'))
        limVal = 200;
        t = Compliance_CladdingGS(M.Nodes,M.Trig,M.Pdis)/limVal;
        n = Area3DFunc(M.Nodes,M.Trig);
        f0s{8} = t./n;
        if derivatives
            dt = Compliance_CladdingGS(M.Nodes,M.Trig,M.Pdis,M.dNodesdx)/limVal;
            dn = Area3DFunc(M.Nodes,M.Trig,M.dNodesdx);
            df0s{8} = (dt.*repmat(n,1,size(dt,2))-repmat(t,1,size(dn,2)).*dn)./repmat(n,1,size(dt,2)).^2;
        end      
    end 
    %_________________________________________________________Eulerbuckling
    if any(strcmp(p.Results.func,'Eulerbuckling'))
        [L,dLdx] = elemsizes(M.Nodes,M.Elements,M.Types,M.dNodesdx);
        Matprop = @(col) M.Materials(M.Elements(:,4),col);
        E = Matprop(2);
        Iy = prop(6);   Iz = prop(7);
        Ncr = @(Iax) pi^2*E.*Iax./(L.^2);
        buckleFunc   = @(Iax) N./Ncr(Iax);
        func = [buckleFunc(Iz);buckleFunc(Iy)];
        [fCluster,dfCluster] = clusterFunc(func,10,M.info.shift);
       
        f0s{11} = fCluster;
        if derivatives
            dIy = dprop(6); dIz = dprop(7);
            dNcr = @(Iax,dIax) pi^2*E.*dFrac(Iax,dIax,L.^2,2*repmat(L,1,nVar).*dLdx,nVar);
            fN = @(N,dN,Nc,dNc) dFrac(N,dN,Nc,dNc,nVar);
            dbuckleFunc = @(Iax,dIax) fN(N, dN, Ncr(Iax),dNcr(Iax,dIax));
            df0s{11} = dfCluster*[dbuckleFunc(Iz,dIz);dbuckleFunc(Iy,dIy)];
        end
    end   
    %_________________________________________________________domain
    if any(strcmp(p.Results.func,'domain'))
        Perimeter = get_perimeter(M.Nodes,M.Elements,M.info.domain);
        NdVrt = cell2mat(Perimeter(:,1));
        NdEdg = cell2mat(Perimeter(:,2));
        NdInt = M.Nodes(setdiff(M.Nodes(:,1),[unique(M.Elements(:,7));NdVrt(:);NdEdg]),:);
        nI    = size(NdInt,1);
        nVrt = size(NdVrt,1);
         Con  = cell(nVrt*nI,1);
        dCon  = cell(nVrt*nI,1);
        for i=1:size(NdVrt,1)
            n1 = NdVrt(i,1); n2 = NdVrt(i,2);
            for ndI = 1:nI
                 A = Area3DFunc(M.Nodes,[1,n1,n2,NdInt(i)]);
                dA = Area3DFunc(M.Nodes,[1,n1,n2,NdInt(i)],M.dNodesdx);
                ind = (i-1)*nI+ndI;
                 Con{ind} = -sign(A)*A+1;
                dCon{ind} = -sign(A)*dA;
            end
        end
        f0s{12} = cell2mat(Con);
        if derivatives
            df0s{12} = cell2mat(dCon);
        end
    end 
    %_______________________________________________________________LoadLim
    if any(strcmp(p.Results.func,'loadLim'))
        Perimeter = Problem.getPeri(M.Nodes,M.Elements,M.info.domain);
        NdVrt = unique(cell2mat(Perimeter(:,1)));
        try 
            Pinit = Problem.Pinit;
        catch
            Problem.Pinit = P(NdVrt,lc);
            Pinit = Problem.Pinit;
        end
        f0s{13} = P(NdVrt,lc)./Pinit;
        if derivatives
            df0s{13} = dP(NdVrt,:,lc)./repmat(Pinit,1,size(dP,2));
        end
    end
    
    %__________________________________________________________________Hmax
    if any(strcmp(p.Results.func,'Hmax'))
        validID = setdiff(M.Nodes(:,1),M.Elements(:,7));
        Hmax = 4;
        f0s{14} = M.Nodes(validID,4)/Hmax;
        if derivatives
            df0s{14} = reshape(M.dNodesdx(validID,4,:),size(validID,1),size(M.dNodesdx,3))/Hmax;
        end
    end
 f0CELL{lc} = cell2mat( f0s);
if derivatives
    df0CELL{lc} = cell2mat(df0s);
end
end

 f = cell2mat( f0CELL);
df = cell2mat(df0CELL);



if strcmp(p.Results.ObjCon,'objective')
    varargout{1}=f(1);
    if derivatives
        if nargin > 8 && varargin{1} && M.info.shift
            key = zVar_key(Problem.SID);
            df(:,key.indX) = regularizator(M.Nodes(key.indC,:),df(:,key.indX),Problem.R);
        end 
        varargout{2}=df(1,:);
    end  
else
    varargout{1}=f-1; varargout{2}=[];
    if derivatives
        if nargin > 8 && varargin{1}
            key = zVar_key(Problem.SID);
            df(:,key.indX) = regularizator(M.Nodes(key.indC,:),df(:,key.indX),Problem.R);
        end 
        varargout{3}=df'; varargout{4}=[];
    end  
end
end

%% Auxiliary functions

function [A] = AreaFunc(Nodes,Trigs,varargin)
nTrig = size(Trigs,1);
if nargin < 3
    nC = 1;
    dNd = 1;
else
    nC = 3*size(Nodes,1);
    dNodesdx = varargin{1};
    dNd = reshape(dNodesdx(:,2:4,:),3*size(Nodes,1),size(dNodesdx,3));
end
Av = zeros(nTrig,nC);
for i = 1:nTrig
    NdT  = Trigs(i,2:4);
    Cor  = Nodes(NdT,2:4);
    x1 = Cor(1); y1 = Cor(4);
    x2 = Cor(2); y2 = Cor(5);
    x3 = Cor(3); y3 = Cor(6);
    if nargin < 3
        j = 1;
        Az = 0.5*((x1-x2)*(y1-y3)-(x1-x3)*(y1-y2));
    else
        j = Addkron(NdT,size(Nodes,1)*(0:2));  
        Az = 0.5*[y2-y3,-x2+x3,0,-y1+y3,x1-x3,0,y1-y2,-x1+x2,0];
    end
    Av(i,j) = Az;
end
[Ai,Aj] = ndgrid((1:nTrig),(1:nC));
A = sparse(Ai(:),Aj(:),Av(:),nTrig,nC)*dNd;
end




function [A] = Area3DFunc(Nodes,Trigs,varargin)
nTrig = size(Trigs,1);
if nargin < 3
    nC = 1;
    dNd = 1;
else
    nC = 3*size(Nodes,1);
    dNodesdx = varargin{1};
    dNd = reshape(dNodesdx(:,2:4,:),3*size(Nodes,1),size(dNodesdx,3));
end
Av = zeros(nTrig,nC);
for i = 1:nTrig
    NdT  = Trigs(i,2:4);
    Cor  = Nodes(NdT,2:4);
    x1 = Cor(1); y1 = Cor(4); z1 = Cor(7);
    x2 = Cor(2); y2 = Cor(5); z2 = Cor(8);
    x3 = Cor(3); y3 = Cor(6); z3 = Cor(9);
    Ax = ((y1-y2)*(z1-z3)-(y1-y3)*(z1-z2));
    Ay = ((z1-z2)*(x1-x3)-(z1-z3)*(x1-x2));
    Az = ((x1-x2)*(y1-y3)-(x1-x3)*(y1-y2));
    Norm = (Ax^2+Ay^2+Az^2)^(0.5);
    if nargin < 3
        j = 1;
        A = 0.5*Norm;
    else
        j = Addkron(NdT,size(Nodes,1)*(0:2));  
        dNorm = [((y2-y3)*Az-(z2-z3)*Ay)/Norm,(-(x2-x3)*Az+(z2-z3)*Ax)/Norm,((x2-x3)*Ay-(y2-y3)*Ax)/Norm,...
                (-(y1-y3)*Az+(z1-z3)*Ay)/Norm,((x1-x3)*Az-(z1-z3)*Ax)/Norm,(-(x1-x3)*Ay+(y1-y3)*Ax)/Norm,...
                 ((y1-y2)*Az-(z1-z2)*Ay) /Norm,(-(x1-x2)*Az+(z1-z2)*Ax)/Norm,((x1-x2)*Ay-(y1-y2)*Ax)/Norm];
        A = 0.5*dNorm;
    end
    Av(i,j) = A;
end
[Ai,Aj] = ndgrid((1:nTrig),(1:nC));
A = sparse(Ai(:),Aj(:),Av(:),nTrig,nC)*dNd;
end




function [S] = SurfaceFunc(Nodes,NID,varargin)
nNd = size(NID,1);
x = Nodes(NID,2); y = Nodes(NID,3); z = Nodes(NID,4);
w = 17;
if nargin < 3
    nC = 1;
    dNd = 1;
    Si = NID;
    Sj = ones(nNd,1);
    Sv = (z-4.*(w^(-4)).*(w-x).*(w+x).*(w-y).*(w+y)).^2;
else
    dNodesdx = varargin{1};
    nC = 3*size(Nodes,1);
    dNd = reshape(dNodesdx(:,2:4,:),3*size(Nodes,1),size(dNodesdx,3));
    Si = repmat(NID,3,1);
    Sj = Addkron(size(Nodes,1).*(0:2)',NID);
    dfdx = 16*x.*(w-y).*(w+y).*(w^4*z-4*(w-x).*(w+x).*(w-y).*(w+y))/w^8;
    dfdy = 16*y.*(w-x).*(w+x).*(w^4*z-4*(w-x).*(w+x).*(w-y).*(w+y))/w^8;
    dfdz = 2*(w^4*z-4*(w-x).*(w+x).*(w-y).*(w+y))/w^4;
    Sv = [dfdx;dfdy;dfdz];
end
S = sparse(Si(:),Sj(:),Sv(:),nNd,nC)*dNd;
end


function [key] = zVar_key(SID)
aux = (1:size(SID,1))';
key = struct;
key.indX = aux(mod(SID,3)==0);
key.indC = SID(mod(SID,3)==0)/3;
end

function [dfW] = regularizator(Nds,df,R)
[x,y] = ndgrid(Nds(:,2),Nds(:,3));
dij = ((x-x').^2+(y-y').^2).^0.5;
Hij = max(R-dij,0);
dfW = (df*Hij)./(ones(size(df))*Hij);
end

function [fCluster,dfCluster] = clusterFunc(func,varargin)
nElm = size(func,1);
p = 10;
if nargin <2
    nGroup = min(10,nElm);
else
    nGroup = varargin{1};
end
Paux = @(x,p) x.^p;
elemPerGroup = floor(nElm/nGroup);
fCluster = zeros(nGroup,1);
dfCluster = zeros(nGroup,nElm);
for i=1:nGroup
    i1 = elemPerGroup*(i-1);
    ind = (i1+1:min(i1+elemPerGroup,nElm));
    group = func(ind,:);
    if sum(group(:))<0 && varargin{2}
        eps = -2*min(group(:));
    else
        eps = 0;
    end
    fCluster(i) = real(sum(Paux(group+eps,p),1).^(1/p))-eps;
    dfCluster(i,ind) = real(sum(Paux(group+eps,p),1).^(1/p-1)*Paux(group+eps,p-1));
end
end



