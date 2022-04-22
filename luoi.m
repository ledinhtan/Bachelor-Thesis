clear all

Nx=2; % So doan tren truc x
Ny=2; % So doan tren truc y
nx=Nx+1; % So node tren truc x
ny=Ny+1; % So node tren truc y

hx=1/Nx;
hy=1/Ny;

% toa do nut
for j=1:ny
    for i=1:nx
        node=(j-1)*nx+i;
        coordinate(node,1)=(i-1)*hx;
        coordinate(node,2)=(j-1)*hy;
    end
end

% luoi
for j=1:ny-1
    count=1;
    for i=1:nx-1
        e=2*(j-1)*(nx-1)+count;
        element(e,1)=(j-1)*nx+i;
        element(e,2)=(j-1)*nx+i+1;
        element(e,3)=j*nx+i+1;
        element(e+1,1)=(j-1)*nx+i;
        element(e+1,2)=j*nx+i+1;
        element(e+1,3)=j*nx+i;
        count=count+2;
    end
end
%coordinate
%element
[nodes2element,nodes2edge,noedges,edge2element,interioredge,exterioredge]=edge(element,coordinate);

