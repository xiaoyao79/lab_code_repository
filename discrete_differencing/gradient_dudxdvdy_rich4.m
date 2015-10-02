function [dudx,dudy,dvdx,dvdy] = gradient_dudxdvdy_rich4(u,v,dx,dy)
% 4th order richardson gradient scheme
% only does the central differences
[row,col] = size(u);

dudx=zeros(row,col);
dudy=zeros(row,col);
dvdx=zeros(row,col);
dvdy=zeros(row,col);
dwdx=zeros(row,col);
dwdy=zeros(row,col);

A=[1239 272 1036 0 -69];
k=[1 2 4 8];
for i=9:row-8
    for j=9:col-8
        cdudx=0;cdvdx=0;cdudy=0;cdvdy=0;cdwdx=0;cdwdy=0;
        for m=1:4;
            cdudx=A(m+1)*(u(i,j+k(m))-u(i,j-k(m)))/(2*k(m)*dx)+cdudx;
            cdudy=A(m+1)*(u(i+k(m),j)-u(i-k(m),j))/(2*k(m)*dy)+cdudy;
            cdvdx=A(m+1)*(v(i,j+k(m))-v(i,j-k(m)))/(2*k(m)*dx)+cdvdx;
            cdvdy=A(m+1)*(v(i+k(m),j)-v(i-k(m),j))/(2*k(m)*dy)+cdvdy;
%             cdwdx=A(m+1)*(w(i+k(m),j)-w(i-k(m),j))/(2*k(m)*dx)+cdwdx;
%             cdwdy=A(m+1)*(w(i+k(m),j)-w(i-k(m),j))/(2*k(m)*dy)+cdwdy;

        end
        dudx(i,j)=1/A(1)*cdudx;
        dudy(i,j)=1/A(1)*cdudy;
        dvdx(i,j)=1/A(1)*cdvdx;
        dvdy(i,j)=1/A(1)*cdvdy;
%         dwdx(i,j)=1/A(1)*cdwdx;
%         dwdy(i,j)=1/A(1)*cdwdy;
        
    end
end
end