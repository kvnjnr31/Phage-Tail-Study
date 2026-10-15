function [xunit,yunit] = circle(x,y,r,cir_pts)

th = linspace(0,2*pi,cir_pts);
xunit = r* cos(th) + x;
yunit = r* sin(th) + y;

