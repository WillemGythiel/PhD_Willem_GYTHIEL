function [Nodes,Elements,dofH,PH] = domeKiewitt(submK,subcK,r)  
q = 2000;               % gravity load [N/m²]
alpha_mK = pi/2/submK;
alpha_cK = pi*2/subcK;
alpha_ctK = alpha_cK/submK;
%% NODES
    % Determining the placement of the nodes of the different dome types
    % Schwedler dome
    
        Nodes = [];
        k = 1;
        for j = 0:submK-1
            Y(j+1) = cos(j*alpha_mK)*r;
            for jj = 0:subcK*submK-1-j*subcK
                z(j+1,jj+1) = sin(j*alpha_mK)*r;
                y(j+1,jj+1) = cos(jj*alpha_ctK)*Y(j+1);
                x(j+1,jj+1) = sin(jj*alpha_ctK)*Y(j+1);

                Nodes = [Nodes; k x(j+1,jj+1) y(j+1,jj+1) z(j+1,jj+1)];
                k = k+1;
            end
            alpha_ctK = pi*2/(subcK*(submK-j-1));
        end
        Nodes = [Nodes; 
                    k 0 0 r;
                    k+1 0 0 0];
        nNodeKiew = Nodes(end,1);
        
        Nodes(:,2) = Nodes(:,2)+20*ones(size(Nodes(:,2)));
                      

    %% ELEMENTS [iElt iType iSec iMat iNode1 iNode2 iNode3]
    % Defining the elements using the previously determined nodes
   
        subcK2 = subcK*submK;

        Elements = [1 1 1 1 1 2 nNodeKiew];                                             % circumferential layers
        Elements = reprow(Elements,1,subcK2-1,[1 0 0 0 1 1 0]);
        Elements(end,6) = Elements(1,5);               
        nElemK = Elements(end,1);
        Elements = reprow(Elements,nElemK-subcK2+1:nElemK-subcK,1,[nElemK 0 0 0 subcK2 subcK2 0]);
        Elements(end,6) = Elements(nElemK+1,5);
        nElemK2 = Elements(end,1);
        for i = 1:submK-2
            Elements = reprow(Elements,nElemK+1:nElemK2-subcK,1,[subcK2-(subcK*i) 0 0 0 subcK2-(subcK*i) subcK2-(subcK*i) 0]);
            Elements(end,6) = Elements(nElemK2+1,5);
            nElemK = nElemK2;
            nElemK2 = Elements(end,1);
        end
        nElemK = nElemK2;
        Elements = [Elements;                                                       % meridionals  
                        nElemK+1 1 1 1 1 1+subcK2 nNodeKiew;
                        nElemK+2 1 1 1 2 1+subcK2 nNodeKiew];
        Elements = reprow(Elements,nElemK+1:nElemK+2,submK-1,[2 0 0 0 1 1 0]);
        nElemK = Elements(end,1);
        Elements = reprow(Elements,nElemK2+2:nElemK,subcK-1,[submK*2-1 0 0 0 submK submK-1 0]);
        Elements(end,:) = []; 
        Elements(end,6) = Elements(1+subcK2,5);
        nElemK3 = Elements(end,1);
        for i = 0:submK-2
            Elements = reprow(Elements,nElemK2+1:nElemK-2,1,[nElemK3-nElemK2 0 0 0 subcK2-(subcK*i) subcK2-(subcK*(i+1)) 0]);
            nElemK = Elements(end,1);
            nElemK2 = nElemK3;
            Elements = reprow(Elements,nElemK2+2:nElemK,subcK-1,[nElemK-(nElemK2+2)+1 0 0 0 submK-1-i submK-2-i 0]);
            Elements(end,:) = []; 
            Elements(end,6) = Elements(end,5)+1;
            nElemK3 = Elements(end,1);
        end
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

  
        % ElementsKiewShell=[EltID TypID SecID MatID n1 n2 n3 n4]
        add1KShell = 0;
        add2KShell = 0;
        ElementsShell = [];
        for j = 1:submK-1
            if isempty(ElementsShell) 
                NElemKiew1 = 0; 
            else
                NElemKiew1 = length(ElementsShell(:,1));
            end
            ElementsShell = [ElementsShell; NElemKiew1+1 1 1 1 add1KShell+2 add1KShell+1 add2KShell+submK*subcK+1 add1KShell+2];
            ElementsShell = reprow(ElementsShell,NElemKiew1+1,submK-j,[1 0 0 0 1 1 1 1]);
            NElemKiew2 = length(ElementsShell(:,1));
            ElementsShell = reprow(ElementsShell,NElemKiew1+1:NElemKiew2,subcK-1,[submK-j+1 0 0 0 submK-j+1 submK-j+1 submK-j submK-j+1]);
            ElementsShell(end,:) = ElementsShell(end,:) - [0 0 0 0 subcK*(submK+1-j) 0 subcK*(submK-j) subcK*(submK+1-j)];
            NElemKiew1 = length(ElementsShell(:,1));
            ElementsShell = [ElementsShell; NElemKiew1+1 1 1 1 add1KShell+2 add2KShell+submK*subcK+1 add2KShell+submK*subcK+2 add1KShell+2];
            ElementsShell = reprow(ElementsShell,NElemKiew1+1,submK-1-j,[1 0 0 0 1 1 1 1]);
            NElemKiew2 = length(ElementsShell(:,1));
            ElementsShell = reprow(ElementsShell,NElemKiew1+1:NElemKiew2,subcK-1,[submK-j 0 0 0 submK-j+1 submK-j submK-j submK-j+1]);
            ElementsShell(end,:) = ElementsShell(end,:) - [0 0 0 0 0 0 subcK*(submK-j) 0];
            add1KShell = add1KShell + (submK-j+1)*subcK;
            add2KShell = add2KShell + (submK-j)*subcK;
        end
        NElemKiew1 = length(ElementsShell(:,1));
        ElementsShell = [ElementsShell; NElemKiew1+1 1 1 1 2+subcK*(((submK+1)*submK)/2-1) 1+subcK*(((submK+1)*submK)/2-1) 1+subcK*(((submK+1)*submK)/2) 2+subcK*(((submK+1)*submK)/2-1)]; 
        ElementsShell = reprow(ElementsShell, NElemKiew1+1, subcK-2, [1 0 0 0 1 1 0 1]);
        NElemKiew1 = length(ElementsShell(:,1));
        ElementsShell = [ElementsShell; NElemKiew1+1 1 1 1 1+subcK*(((submK+1)*submK)/2-1) subcK*(((submK+1)*submK)/2) 1+subcK*(((submK+1)*submK)/2) 1+subcK*(((submK+1)*submK)/2-1)];

        ElemLoad = accel_shell4(Accelxyz,Nodes,ElementsShell,SectionsShell,MaterialsShell);
        PI = elemloads(ElemLoad,Nodes,ElementsShell,TypesShell,dofI);
        PH = elemloads(ElemLoad,Nodes,ElementsShell,TypesShell,dofH);

end

  