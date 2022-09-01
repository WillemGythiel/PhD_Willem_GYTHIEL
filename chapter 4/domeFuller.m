function [Nodes,Elements,dofH,PH] = domeFuller(subF,r)  
q = 2000;               % gravity load [N/m²]
    kF = 1;
    rotkF = 1;
%% NODES
        Nodes = [];
        % 5 vertices of the icosahedron
        VerticesIco = [                             0,                                      0,                     r;  % v1
                                r*sin(pi/2-atan(1/2)),                                      0, r*cos(pi/2-atan(1/2));  % v2
                    r*cos(8*pi/5)*sin(pi/2-atan(1/2)),      r*sin(8*pi/5)*sin(pi/2-atan(1/2)), r*cos(pi/2-atan(1/2));  % v3
                    r*cos(9*pi/5)*sin(pi/2+atan(1/2)),      r*sin(9*pi/5)*sin(pi/2+atan(1/2)), r*cos(pi/2+atan(1/2));  % v4
                    r*cos(7*pi/5)*sin(pi/2+atan(1/2)),      r*sin(7*pi/5)*sin(pi/2+atan(1/2)), r*cos(pi/2+atan(1/2))]; % v5
        % vertices on 5 edges           
        for j = 0:subF-1
            VerticesEdgev2v1(j+1,:) = VerticesIco(2,:)+j/subF*(VerticesIco(1,:)-VerticesIco(2,:));
            VerticesEdgev3v1(j+1,:) = VerticesIco(3,:)+j/subF*(VerticesIco(1,:)-VerticesIco(3,:));
            VerticesEdgev4v2(j+1,:) = VerticesIco(4,:)+j/subF*(VerticesIco(2,:)-VerticesIco(4,:));
            VerticesEdgev4v3(j+1,:) = VerticesIco(4,:)+j/subF*(VerticesIco(3,:)-VerticesIco(4,:));
            VerticesEdgev5v3(j+1,:) = VerticesIco(5,:)+j/subF*(VerticesIco(3,:)-VerticesIco(5,:));
        end
        % all vertices on icosahedron
        for j = subF/2:subF-1
            for jj = 0:j-1     % Vertices triangle v2v3v4         
                Nodes = [Nodes; VerticesEdgev4v2(j+1,:)+jj/j*(VerticesEdgev4v3(j+1,:)-VerticesEdgev4v2(j+1,:))];
                kF = kF+1;
            end
            for jj = 0:subF-1-j % Vertices triangle v3v4v5
                Nodes = [Nodes; VerticesEdgev4v3(j+1,:)+jj/(subF-j)*(VerticesEdgev5v3(j+1,:)-VerticesEdgev4v3(j+1,:))];
                kF = kF+1;
            end
            for rot = 1:4      % Rotation of vertices triangle v2v3v4 and v3v4v5 
                Nodes = [Nodes; Nodes(rotkF:kF-1,:)*[cos(rot*2*pi/5), -sin(rot*2*pi/5), 0; sin(rot*2*pi/5), cos(rot*2*pi/5), 0; 0, 0, 1]];
            end
            kF = kF+4*(kF-rotkF);
            rotkF = kF;
        end
        for j = 0:subF-1
            for jj = 0:subF-1-j  % Vertices triangle v1v2v3
                Nodes = [Nodes; VerticesEdgev2v1(j+1,:)+jj/(subF-j)*(VerticesEdgev3v1(j+1,:)-VerticesEdgev2v1(j+1,:))];
                kF = kF+1;
            end
            for rot = 1:4       % Rotation of vertices triangle v1v2v3
                Nodes = [Nodes; Nodes(rotkF:kF-1,:)*[cos(rot*2*pi/5), -sin(rot*2*pi/5), 0; sin(rot*2*pi/5), cos(rot*2*pi/5), 0; 0, 0, 1]];
            end
            kF = kF+4*(kF-rotkF);
            rotkF = kF;
        end
        Nodes = [Nodes; VerticesIco(1,:)]; % addition of top node
        % projection on sphere
        for j = 1:kF
            Nodes(j,:) = Nodes(j,:)*r/norm(Nodes(j,:));
        end 
        Nodesk = [1:kF];
        Nodes = [Nodesk.', Nodes];
        kF = kF+1;
        Nodes = [Nodes; kF, 0, 0, 0];

        Nodes(:,2) = Nodes(:,2)-20*ones(size(Nodes(:,2)));                   
    %% ELEMENTS [iElt iType iSec iMat iNode1 iNode2 iNode3]
    % Defining the elements using the previously determined nodes
    

        Elements = [1 1 1 1 1 2 kF];
        Elements = reprow(Elements, 1, 5*subF-2, [1 0 0 0 1 1 0]);
        Elements = [Elements; 5*subF 1 1 1 5*subF 1 kF];
        Elements = reprow(Elements, 1:5*subF, subF/2, [5*subF 0 0 0 5*subF 5*subF 0]);
        add = 0;
        for j = 1:subF-1
            Elements = [Elements; add+5*subF*(subF/2+1)+1 1 1 1 add+5*subF*(subF/2+1)+1 add+5*subF*(subF/2+1)+2 kF];   
            Elements = reprow(Elements, add+5*subF*(subF/2+1)+1, 5*(subF-j)-2, [1 0 0 0 1 1 0]);
            addNew = add + (subF-j)*5;
            Elements = [Elements; addNew+5*subF*(subF/2+1) 1 1 1 addNew+5*subF*(subF/2+1) add+5*subF*(subF/2+1)+1 kF];
            add = addNew;
        end
        NElem1 = length(Elements);
        Elements = [Elements; NElem1+1 1 1 1 1 5*subF+1 kF];
        Elements = reprow(Elements, NElem1+1, 5*subF*subF/2-1, [1 0 0 0 1 1 0]);
        NElem1 = length(Elements);
        Elements = [Elements; NElem1+1 1 1 1 1 5*subF+2 kF];
        Elements = reprow(Elements, NElem1+1, 5*subF-2, [1 0 0 0 1 1 0]);
        NElem2 = length(Elements);
        Elements = [Elements; NElem2+1 1 1 1 5*subF 5*subF+1 kF];
        Elements = reprow(Elements, NElem1+1:NElem2+1, subF/2-1, [5*subF 0 0 0 5*subF 5*subF 0]);
        add1 = 0;
        add2 = 0;
        for j = 1:subF-1
            NElem2 = length(Elements);
            Elements = [Elements; NElem2+1 1 1 1 add1+1+subF*5*subF/2 add2+1+subF*5*(subF/2+1) kF;
                                  NElem2+2 1 1 1 add1+2+subF*5*subF/2 add2+1+subF*5*(subF/2+1) kF];
            Elements = reprow(Elements,NElem2+1:NElem2+2,subF-j,[2 0 0 0 1 1 0]);
            NElem1 = length(Elements);
            Elements = reprow(Elements,NElem2+2:NElem1,5-1,[2*(subF-j)+1 0 0 0 subF-j+1 subF-j 0]);
            Elements(end,:) = []; 
            Elements(end,6) = Elements(add2+1+subF*5*(subF/2+1),5);
            add1 = add1 + (subF-j+1)*5;
            add2 = add2 + (subF-j)*5;
        end
        NElem2 = length(Elements);
        Elements = [Elements; NElem2+1 1 1 1 add1+1+subF*5*subF/2 kF-1 kF];
        Elements = reprow(Elements,NElem2+1,5-1,[1 0 0 0 1 0 0]);
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

    
        % ElementsFullShell=[   EltID TypID SecID MatID     n1      n2      n3    n4]
        ElementsShell = [       1     1     1     1      1 subF*5+1 subF*5+2     1;
                                    2     1     1     1      2       1 subF*5+2     2];
        ElementsShell = reprow(ElementsShell,      1:2, subF*5-2, [     2 0 0 0     1     1     1     1]);
        ElementsShell = [                            ElementsShell;
                             subF*10-1     1     1     1  subF*5  subF*10 subF*5+1 subF*5;
                               subF*10     1     1     1      1   subF*5 subF*5+1     1];
        ElementsShell = reprow(ElementsShell, 1:subF*10, subF/2-1, [subF*10 0 0 0 subF*5 subF*5 subF*5 subF*5]); 
        add1FShell = 0;
        add2FShell = 0;
        subcF = 5;
        for j = 1:subF-1
            NElemFullShell1 = length(ElementsShell(:,1));
            ElementsShell = [ElementsShell; NElemFullShell1+1 1 1 1 add1FShell+2+5*subF*subF/2 add1FShell+1+5*subF*subF/2 add2FShell+1+5*subF*(subF/2+1) add1FShell+2+5*subF*subF/2];
            ElementsShell = reprow(ElementsShell,NElemFullShell1+1,subF-j,[1 0 0 0 1 1 1 1]);
            NElemFullShell2 = length(ElementsShell(:,1));
            ElementsShell = reprow(ElementsShell,NElemFullShell1+1:NElemFullShell2,subcF-1,[subF-j+1 0 0 0 subF-j+1 subF-j+1 subF-j subF-j+1]);
            ElementsShell(end,:) = ElementsShell(end,:) - [0 0 0 0 subcF*(subF+1-j) 0 subcF*(subF-j) subcF*(subF+1-j)];
            NElemFullShell1 = length(ElementsShell(:,1));
            ElementsShell = [ElementsShell; NElemFullShell1+1 1 1 1 add1FShell+2+5*subF*subF/2 add2FShell+1+5*subF*(subF/2+1) add2FShell+2+5*subF*(subF/2+1) add1FShell+2+5*subF*subF/2];
            ElementsShell = reprow(ElementsShell,NElemFullShell1+1,subF-1-j,[1 0 0 0 1 1 1 1]);
            NElemFullShell2 = length(ElementsShell(:,1));
            ElementsShell = reprow(ElementsShell,NElemFullShell1+1:NElemFullShell2,subcF-1,[subF-j 0 0 0 subF-j+1 subF-j subF-j subF-j+1]);
            ElementsShell(end,:) = ElementsShell(end,:) - [0 0 0 0 0 0 subcF*(subF-j) 0];
            add1FShell = add1FShell + (subF-j+1)*subcF;
            add2FShell = add2FShell + (subF-j)*subcF;
        end
        NElemFullShell1 = length(ElementsShell(:,1));
        ElementsShell = [ElementsShell; NElemFullShell1+1 1 1 1 2+kF-7 1+kF-7 6+kF-7 2+kF-7]; 
        ElementsShell = reprow(ElementsShell, NElemFullShell1+1, subcF-2, [1 0 0 0 1 1 0 1]);
        NElemFullShell1 = length(ElementsShell(:,1));
        ElementsShell = [ElementsShell; NElemFullShell1+1 1 1 1 1+kF-7 5+kF-7 6+kF-7 1+kF-7];

        ElemLoad = accel_shell4(Accelxyz,Nodes,ElementsShell,SectionsShell,MaterialsShell);
        PI = elemloads(ElemLoad,Nodes,ElementsShell,TypesShell,dofI);
        PH = elemloads(ElemLoad,Nodes,ElementsShell,TypesShell,dofH);
        
end