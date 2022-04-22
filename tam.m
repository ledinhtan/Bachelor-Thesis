%for i=1:sdkb
%   B(i,:)=0.0; K(i,:)=0.0; M(i,:)=0.0;
%   B(:,i)=0.0; K(:,i)=0.0; M(:,i)=0.0;
   
%   B(i,i)=1.0; K(i,i)=1.0; M(i,i)=1.0;

%   B(i,:)=[]; K(i,:)=[]; M(i,:)=[]; 
%   B(:,i)=[]; K(:,i)=[]; M(:,i)=[];
%end
%A=[-beta*B, -(K+alpha*B); M, zeros(noedges-sdkb)];
%D=[M,zeros(noedges-sdkb); zeros(noedges-sdkb),M];


A=[-beta*B, -(K+alpha*B); M, zeros(noedges)];
D=[M,zeros(noedges); zeros(noedges),M];

for i=1:sdkb
   A(i,:)=0.0; A(:,i)=0.0; A(i,i)=1.0;  
   A(i+noedges,:)=0.0; A(:,i+noedges)=0.0; A(i+noedges,i+noedges)=1.0;  
   D(i,:)=0.0; D(:,i)=0.0; D(i,i)=1.0;  
   D(i+noedges,:)=0.0; D(:,i+noedges)=0.0; D(i+noedges,i+noedges)=1.0;  
end

