function [dKeG] = DiffKeGCS_Frame_Coor(Coor,SecProp,MatProp,d)
 Te  = Te_Frame(Coor);
dTe  = Te_Frame(Coor,'shape',d);
 KeL = KeL_Frame(Coor,SecProp,MatProp);
dKeL = KeL_Frame(Coor,SecProp,MatProp,'shape',d);
dKeG = dTe'*KeL*Te+Te'*dKeL*Te+Te'*KeL*dTe;
end  
