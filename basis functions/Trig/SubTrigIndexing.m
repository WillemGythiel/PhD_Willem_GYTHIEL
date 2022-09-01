function [TrigsNodeNr] = SubTrigIndexing(n)
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



% NupTrigs = (n+1)*n/2;
% NloTrigs = (n-1)*n/2;
% % Indexing of the upper triangles
% UpTrigs = zeros(NupTrigs,3);
% for i = 1:n
%     Vertex1 = (i*(i-1)/2+1:i*(i+1)/2)';
%     Index   = Vertex1;
%     UpTrigs(Index,:) = [Vertex1, Vertex1+i, Vertex1+i+1];
% end
% % Indexing of the lower triangles
% LoTrigs = zeros(NloTrigs,3);
% if n > 1
% for j = 2:n
%     Vertex1 = ((j-1)*j/2+1:(j-1)*(j+2)/2)';
%     Index   = ((j-2)*(j-1)/2+1:(j-1)*j/2)';
%     LoTrigs(Index,:) = [Vertex1+j+1, Vertex1+1, Vertex1];
% end
% end
% TrigsNodeNr = [UpTrigs; LoTrigs];
