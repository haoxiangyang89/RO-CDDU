clear;
clc;
m=500;

xCan = linspace(0,5,floor(m));
yCan = linspace(0,5,m);

V = zeros(m,m);
Z1 = zeros(m,m);
Z2 = zeros(m,m);
for i = 1:m
    for j = 1:m
        x = xCan(i);
        y = yCan(j);
        %if ((-6*x+8*y <= 3) && (3*x - y <=3))
            Z1(i,j) = max( [-4*x+2*y +0.5*x*y+6, 2*x-y+0.4*x*y+10,-2*x -2*y+0.1*x*y+20,x+2*y+0.01*x*y+5]);
            %Z2(i,j) = max([ -2*x -y+0.1*x*y+18,-x-2*y+0.15*x*y+15]);
        %end
    end
end
V = Z1+Z2;
figure(1)
 mesh(yCan,xCan,V)
 xlabel('y');
 ylabel('x');
 zlabel('Obj');
 figure(2)
  mesh(yCan,xCan,Z1)
 xlabel('y');
 ylabel('x');
 zlabel('Z1');
  figure(3)
  mesh(yCan,xCan,Z2)
 xlabel('y');
 ylabel('x');
 zlabel('Z2');