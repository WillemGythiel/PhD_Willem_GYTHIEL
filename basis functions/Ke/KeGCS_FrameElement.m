function [Ke] = KeGCS_FrameElement(Coor,SecProp,MatProp)
Te  = Te_Frame(Coor);
KeL = KeL_Frame(Coor,SecProp,MatProp);
Ke  = Te'*KeL*Te;
end
