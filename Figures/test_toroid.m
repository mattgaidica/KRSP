%% // define initial tore
% see https://stackoverflow.com/questions/36157881/drawing-circles-on-a-torus-in-matlab
% close all
ff(1000,800);
am = 1.; rm = 3.;
t = linspace(-pi,pi,1440);
p = linspace(0,2.*pi,366);
[T,P] = meshgrid(p,t);
X = (rm + am.*cos(P)).*cos(T);
Y = (rm + am.*cos(P)).*sin(T);
Z = am.*sin(P);
filtStd_mod = flip(filtStd);
filtStd_mod(1,:) = 0;
filtStd_mod(:,720) = NaN;
CO = ind2rgb(uint16(normalize(filtStd_mod','range')*65535),magma(65535));
hsurf = surf(X,Y,Z,CO,'FaceColor','interp','FaceAlpha',0.85,'EdgeColor','none');
axis equal
grid off
box off
view(-130,43);

RGB = ind2rgb(filtStd',magma);