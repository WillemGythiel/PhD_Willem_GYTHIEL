function [info] = TrigInfoGS(n)
iD = struct; iK = struct; iL = struct; 
%____________________________________________________________Number_of_nodes
     nDOF = 6;
      nNd = (n+1)*(n+2)/2;
     iD.n = nDOF*nNd;
%_____________________________________________________________Special_nodes
        k = (1:n)';
     ind1 = k.*(k-1)/2+1;
     ind2 = k+n*(n+1)/2;
     ind3 = flipud(ind1+2*k);
     EdNd = [ind1;ind2;ind3];
     VrNd = [1 nNd-n nNd];
%______________________________________________________________________DOFS
   AllDOF = (1:iD.n);
   ConDOF = Addkron((1:3),nDOF*(VrNd-1));
   iD.dof = setdiff(AllDOF,ConDOF);
%__________________________________________________Stiffness_matrix_Indices
     ElNd = SubTrigNodeNrs(n);  
[DKTi,DKTj] = ndgrid([3;5;4],[3;5;4]);
[TMDi,TMDj] = ndgrid([1;2;6],[1;2;6]);
     iK.i = Addkron(ElNd,[DKTi(:);TMDi(:)]);
     iK.j = Addkron(ElNd,[DKTj(:);TMDj(:)]);
%______________________________________________________________Load_Indices
   MatMat = repmat((0:n)/n,n+1,1)+1;
       L1 = nonzeros(triu(-MatMat))+2;
       L3 = nonzeros(triu(MatMat'))-1;
       L2 = 1-L1-L3;
     iL.l = [L1 L2 L3];
  TriNdNr = (1:nNd);
     iL.i = Addkron((1:nDOF),6*(TriNdNr-1))';
     iL.j = ones(nDOF*nNd,1);
%___________________________________________________________Load_Multiplier
 LM       = 2*ones(nNd,1);
 LM(EdNd) = 1;
 LM(VrNd) = 1/3;
     iL.f = diag(LM)/n^2;
     iL.v = zeros(nNd,6);
info = struct;
info.iD = iD; info.iK = iK; info.iL = iL; info.n = n;
end

function [TrigsNodeNr] = SubTrigNodeNrs(n)
     AuxMat = repmat((1:n),n,1);
% Indexing of the upper triangles
        Nd1 = (1:n*(n+1)/2)';
        inc = nonzeros(triu(AuxMat));
    UpTrigs = [Nd1, Nd1+inc, Nd1+inc+1];
% Indexing of the lower triangles
        Nd2 = UpTrigs(1:n*(n-1)/2,2);
        add = nonzeros(triu(AuxMat,1));
    LoTrigs = [Nd2+add+1, Nd2+1, Nd2];
TrigsNodeNr = [UpTrigs; LoTrigs];
end