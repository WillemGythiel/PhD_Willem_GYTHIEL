function [out] = SubTriginfoGS(n)
iD = struct; iK = struct; iS = struct; iL = struct; 
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
   ConDOF = [Addkron((1:3),nDOF*(VrNd-1)),3*(EdNd'-1)+6];
   iD.dof = setdiff(AllDOF,ConDOF);
%____________________________________________________________Spring_Indices
     DOFk = [(1:6)',1e5*ones(6,1)];                                            % spring orientation and k value
     iS.i = Addkron(nDOF*(EdNd-1),DOFk(:,1));
     iS.j = iS.i;
     iS.v = repmat(DOFk(:,2),size(iS.i,1)/size(DOFk,1),1);  
%___________________________________________________Stiffness_matrix_Indices
     ElNd = SubTrigIndexing(n);  
      DpE = ones(1,size(ElNd,2)*nDOF);
      pos = Addkron(nDOF*(ElNd-1),(1:nDOF));
       Mi = kron(DpE,pos)';
       Mj = kron(pos,DpE)';
     iK.i = Mi(:);
     iK.j = Mj(:);
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
out = {iD,iK,iS,iL,n};
end
%% AUXILIARY FUNCTIONS
function [TrigsNodeNr] = SubTrigIndexing(n)
     AuxMat = repmat((1:n),n,1);
% Indexing of the upper triangles
        Nd1 = (1:n*(n+1)/2)';
        inc = nonzeros(triu(AuxMat));
    UpTrigs = [Nd1+inc, Nd1+inc+1, Nd1];
% Indexing of the lower triangles
        Nd2 = UpTrigs(1:n*(n-1)/2,2);
        add = nonzeros(triu(AuxMat,1));
    LoTrigs = [Nd2+1, Nd2, Nd2+add+1];
TrigsNodeNr = [UpTrigs; LoTrigs];
end