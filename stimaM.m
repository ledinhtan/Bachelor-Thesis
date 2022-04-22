function M=stimaM(coord,rho)
N=coord(:)*ones(1,3)-repmat(coord,3,1);
D=diag([norm(N([5,6],2)) norm(N([1,2],3)) norm(N([1,2],2))]);
C=spdiags([ones(6,1),ones(6,1),2*ones(6,1),ones(6,1),ones(6,1)],[-4,-2,0,2,4],6,6);
%L=rho*eye(6,6);
%M = D*N'*C*L*N*D/(24*det([1,1,1;coord]));
M = rho*D*N'*C*N*D/(24*det([1,1,1;coord]));