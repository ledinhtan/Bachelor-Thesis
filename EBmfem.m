% EBmfem.m
% The main program
% Update: 24/12/2016
% Mo ta bai toan:
% Tim cac tri rieng ung voi cac mode dao dong cua mot khoang cung hinh chu
% nhat 2D voi mot buc tuong hap thu chua day khong khi. Cac tham so hinh
% hoc cua bai toan la: Chieu dai a = 1.00 m, chieu rong b = 0.75 m. Cac tham
% so vat ly cua bai toan duoc cho nhu sau: 
% Mat do khoi rho = 1 kg/m^3.
% Van toc am thanh c = 340 m/s.
% Cac he so alpha = 5e4 N/m^3 va beta = 200 Ns/m^3.
% Bien - mo ta
% nx ................ Tong so node tren truc x
% ny ................ Tong so node tren truc y
% nn ................ Tong so node 
% Nx ................ So doan tren truc x
% Ny ................ So doan tren truc y
% hx ................ Do dai mot doan tren truc x
% hy ................ Do dai mot doan tren truc y
% N ................. So lan lam min luoi (refine mesh)
% B ................. Ma tran do cung toan cuc
% K ................. Ma tran do cung toan cuc
% M ................. Ma tran do cung toan cuc
% dcond ............. Dieu kien bien cung
% cond .............. Dieu kien bien hap thu
% signum ............ Ma tran luu dau the hien huong cua canh
% coordinate ........ Toa do node toan cuc
% element ........... Ma tran dinh vi nut phan tu
% nodes2element ..... Ma tran dinh vi so phan tu cua hai node
% nodes2edge ........ Ma tran dinh vi so canh cua hai node
% noedges ........... So canh cua luoi tam giac
% edge2element ...... Ma tran dinh vi phan tu cua canh 
%%
close all
clear all
clc
%% TIEN XU LY
%----------------------------------------------------------------------------------
% Tham so vat ly
%----------------------------------------------------------------------------------
rho=1;
c=340;
alpha=5e4;
beta=200;
%----------------------------------------------------------------------------------
% Tham so hinh hoc
%----------------------------------------------------------------------------------
ax=0.0; 
bx=1.0;
ay=0.0;
by=0.75;
%----------------------------------------------------------------------------------
% Luoi PTHH
%----------------------------------------------------------------------------------
N=2; % So lan refine mesh
Nx=8*N; % So doan tren truc x
Ny=6*N; % So doan tren truc y
nx=Nx+1; % So node tren truc x
ny=Ny+1; % So node tren truc y
nn=nx*ny; % Tong so node 
hx=(bx-ax)/Nx;
hy=(by-ay)/Ny;
%----------------------------------------------------------------------------------
% Tao luoi tu dong
%----------------------------------------------------------------------------------
% Toa do nut
for j=1:ny
    for i=1:nx
        node=(j-1)*nx+i;
        coordinate(node,1)=(i-1)*hx;
        coordinate(node,2)=(j-1)*hy;
    end
end

% Phan tu
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
 
[nodes2element,nodes2edge,noedges,edge2element]=edge(element,coordinate);
%----------------------------------------------------------------------------------
% Dieu kien bien Dirichlet
%----------------------------------------------------------------------------------
% Canh day
nodes=1:nx;
for i=1:Nx
   dcond(i)=nodes2edge(nodes(i),nodes(i+1));
end 

% Canh trai
nodes=1:nx:(ny-1)*nx+1;
for i=1:Ny
   dcond(i+Nx)=nodes2edge(nodes(i),nodes(i+1));
end 

% Canh phai
nodes=nx:nx:ny*nx;
for i=1:Ny
   dcond(i+Nx+Ny)=nodes2edge(nodes(i),nodes(i+1));
end 

sdkb=length(dcond); % So dieu kien bien
%----------------------------------------------------------------------------------
% Dieu kien bien Neumann
%----------------------------------------------------------------------------------
nodes=(ny-1)*nx+(1:nx); 
for i=1:length(nodes)-1
   cond(i,:)=[nodes2edge(nodes(i),nodes(i+1))  coordinate(nodes(i+1),1)-coordinate(nodes(i),1)];
end
%% XU LY
%----------------------------------------------------------------------------------
% Assemble matrices M,K
%----------------------------------------------------------------------------------
M=sparse(noedges, noedges); 
K=sparse(noedges,noedges);
for j=1:size(element,1) % lap theo phan tu j
   coord=coordinate(element(j,:),:)';
   I=diag(nodes2edge(element(j,[2 3 1]),element(j,[3 1 2]))); 
   signum=ones(1,3);     
   signum(find(j==edge2element(I,4)))=-1;
   
   M(I,I)= M(I,I)+diag(signum)*stimaM(coord,rho)*diag(signum); 
   K(I,I)= K(I,I)+diag(signum)*stimaK(coord,rho,c)*diag(signum);

end
%----------------------------------------------------------------------------------
% Assemble matrices B
%----------------------------------------------------------------------------------
B=sparse(noedges,noedges);
for i=1:size(cond,1)
   B(cond(i,1),cond(i,1))=B(cond(i,1),cond(i,1))+cond(i,2);
end
%----------------------------------------------------------------------------------  
% Khu dieu kien Dirichlet
%----------------------------------------------------------------------------------
B1=-beta*B;
C=-(K+alpha*B);
M1=M;
for i=1:sdkb
   B1(i,:)=0.0; B1(:,i)=0.0; %B1(i,i)=1.0;
   C(i,:)=0.0; C(:,i)=0.0;
   M(i,:)=0.0; M(:,i)=0.0;
   M1(i,:)=0.0; M1(:,i)=0.0; M1(i,i)=1.0;
end
A=[B1, C; M,zeros(noedges,noedges)];
D=[M1, zeros(noedges,noedges); zeros(noedges,noedges), M1];
plotmesh
opts.tol = 1e-3; 
eigs(A,D,10,'sr', opts)
% eigs(A,D,14,'sr')


