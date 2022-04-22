K=sparse(noedges,noedges);
for e=1:size(element,1)
   I=diag(nodes2edge(element(e,[2 3 1]),element(e,[3 1 2])));
   signum=ones(1,3);  
   
   signum(find(j==edge2element(I,4)))=-1;
   K(I,I)= K(I,I)+diag(signum)*stimaK(coord,rho,c)*diag(signum);
end
