function shapeF = shapeFuncQuad(nStep)
[xi,et] = meshgrid(2*(0:nStep)/nStep-1,2*(0:nStep)/nStep-1);
xi = xi(:); et = et(:);
N1 = (1-xi).*(1-et)/4;
N2 = (1+xi).*(1-et)/4;
N3 = (1+xi).*(1+et)/4;
N4 = (1-xi).*(1+et)/4;
shapeF = [N1,N2,N3,N4];
end