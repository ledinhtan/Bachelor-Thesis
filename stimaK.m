function K=stimaK(coord,rho,c)
N=coord(:)*ones(1,3)-repmat(coord,3,1);
D=diag([norm(N([5,6],2)) norm(N([1,2],3)) norm(N([1,2],2))]);
       
%L=rho*c^2*eye(3);       
%K = D*L*D/(det([1,1,1;coord]));
K = 2*rho*c^2*D*D/(det([1,1,1;coord]));