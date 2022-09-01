function [Nodes,Elements,dofH,PH] = domeSchwedler(submS,subcS,r)  
q = 2000;               % gravity load [N/m²]
alpha_mS = pi/2/submS;
alpha_cS = pi*2/subcS;
%% NODES
    % Determining the placement of the nodes of the different dome types
  
        Nodes = [];
        k = 1;
        for j = 0:submS-1
            Y(j+1) = cos(j*alpha_mS)*r;
                for jj = 0:subcS-1
                    z(j+1,jj+1) = sin(j*alpha_mS)*r;
                    y(j+1,jj+1) = cos(jj*alpha_cS)*Y(j+1);
                    x(j+1,jj+1) = sin(jj*alpha_cS)*Y(j+1);

                    Nodes = [Nodes; k x(j+1,jj+1) y(j+1,jj+1) z(j+1,jj+1)];
                    k = k+1;
                end
        end
        Nodes = [Nodes; 
                    k 0 0 r;
                    k+1 0 0 0];
        nNodeSchw = Nodes(end,1);

    %% ELEMENTS [iElt iType iSec iMat iNode1 iNode2 iNode3]
    % Defining the elements using the previously determined nodes

        Elements = [1 1 1 1 1 2 nNodeSchw];                                             % circumferential layers
        Elements = reprow(Elements,1,subcS-1,[1 0 0 0 1 1 0]);
        Elements(end,6) = Elements(1,5);               
        nElemS = Elements(end,1);
        Elements = reprow(Elements,1:nElemS,submS-1,[nElemS 0 0 0 subcS subcS 0]);  
        nElemS = Elements(end,1);
        Elements = [Elements;                                                       % meridionals
                        nElemS+1 1 1 1 1 subcS+1 nNodeSchw];                                
        Elements = reprow(Elements,nElemS+1,submS-1,[1 0 0 0 subcS subcS 0]);
        nElemS2 = Elements(end,1);
        Elements = reprow(Elements,nElemS+1:nElemS2,subcS-1,[nElemS2-nElemS 0 0 0 1 1 0]); 
        nElemS = Elements(end,1);
        for i = 1:nElemS
            if Elements(i,6) > nNodeSchw-1
               Elements(i,6) = nNodeSchw-1;
            end
        end
        Elements = [Elements;                                                       % diagonals
                        nElemS+1 1 1 1 1 subcS+2 nNodeSchw]; 
        Elements = reprow(Elements,nElemS+1,submS-2,[1 0 0 0 subcS subcS 0]);
        nElemS2 = Elements(end,1);
        Elements = reprow(Elements,nElemS+1:nElemS2,subcS-2,[nElemS2-nElemS 0 0 0 1 1 0]); 
        nElemS = Elements(end,1);
        Elements = [Elements;                                                       % diagonals
                        nElemS+1 1 1 1 subcS subcS+1 nNodeSchw]; 
        Elements = reprow(Elements,nElemS+1,submS-2,[1 0 0 0 subcS subcS 0]);


        Elements(:,3)=Elements(:,1);

   

    %% DOFS AND BOUNDARY CONDITIONS
    % Modelling the boundary conditions by removing degrees of freedom
        Types = {1 'truss'};
        dofs = getdof(Elements,Types);
        EdgeNodes = Nodes(Nodes(:,4)<0.01,:);
        dofI = removedof(dofs,EdgeNodes(1:end-1,1)+0.03);
        dofH = removedof(dofs,Addkron(EdgeNodes(1:end-1,1),(1:3)'*0.01));
        seldof = [1.01;
                1.02;
                EdgeNodes(round(end/2),1)+0.01];
        dofI = removedof(dofI,seldof);
        nDOF = length(dofH(:,1));
    %% LOAD
    % Calculating the load by defining shell elements on which a gravity load is applied. This load is then converted to nodal loads 
    TypesShell = {1 'shell4'};
    SectionsShell = [1 1];             % h = 1m
    MaterialsShell = [1 1 1 1];        % rho = 1kg/m^3
    Accelxyz = [0 0 q];                % q = F/A = m*a/A = rho*V*a/A = rho*h*A*a/A = a (with rho = 1kg/m^3 and h = 1m)
        % ElementsSchwShell=[                                                            EltID TypID SecID MatID                n1                n2            n3                n4]
        ElementsShell = [                                                                1     1     1     1                 2                 1       subcS+2                 2];
        ElementsShell = reprow(ElementsShell,         1, subcS-2, [                  1     0     0     0                 1                 1             1                 1]);
        ElementsShell = [ElementsShell;                                          subcS     1     1     1                 1             subcS       subcS+1                 1];
        ElementsShell = [ElementsShell;                                        subcS+1     1     1     1                 1           subcS+1       subcS+2                 1];
        ElementsShell = reprow(ElementsShell,   subcS+1, subcS-2, [                  1     0     0     0                 1                 1             1                 1]);
        ElementsShell = [ElementsShell;                                        2*subcS     1     1     1             subcS           2*subcS       subcS+1             subcS];
        ElementsShell = reprow(ElementsShell, 1:2*subcS, submS-2, [            2*subcS     0     0     0             subcS             subcS         subcS             subcS]);
        ElementsShell = [ElementsShell;                            2*subcS*(submS-1)+1     1     1     1 subcS*(submS-1)+2 subcS*(submS-1)+1 subcS*submS+1 subcS*(submS-1)+2];
        ElementsShell = reprow(ElementsShell, 2*subcS*(submS-1)+1, subcS-2, [        1     0     0     0                 1                 1             0                 1]);
        ElementsShell = [ElementsShell;                            2*subcS*(submS-1/2)     1     1     1 subcS*(submS-1)+1       subcS*submS subcS*submS+1 subcS*(submS-1)+1];

        ElemLoad = accel_shell4(Accelxyz,Nodes,ElementsShell,SectionsShell,MaterialsShell);
        PI = elemloads(ElemLoad,Nodes,ElementsShell,TypesShell,dofI);
        PH = elemloads(ElemLoad,Nodes,ElementsShell,TypesShell,dofH);

end
    