function [U0,V0] = vfm(U,V,dx,dy,dvdy)
[s1,s2] = size(U);

U0l = zeros(s1,s2-1);
for i = 1:s2-1
    if i == 1
        U0l(:,i) = 0.5*(dvdy(:,2)+dvdy(:,1))*dx;
    elseif i == 2
        U0l(:,i) = 0.5*(dvdy(:,3)+2*dvdy(:,2)+dvdy(:,1))*dx;
    else
        U0l(:,i) = 0.5*(dvdy(:,i+1)+2*sum(dvdy(:,2:i),2)+dvdy(:,1))*dx;
    end
end

U0l= U0l-repmat(U(:,1),[1 size(U0l,2)]);
U0l= cat(2,-U(:,1),U0l);

U0r  = zeros(s1,s2-1);
dvdyflip = fliplr(-dvdy);
for i = 1:s2-1
    if i == 1
        U0r(:,i) = -0.5*(dvdyflip(:,2)+dvdyflip(:,1))*dx;
    elseif i == 2
        U0r(:,i) = -0.5*(dvdyflip(:,3)+2*dvdyflip(:,2)+dvdyflip(:,1))*dx;
    else
        U0r(:,i) = -0.5*(dvdyflip(:,i+1)+2*sum(dvdyflip(:,2:i),2)+dvdyflip(:,1))*dx;
    end
end

U0r= U0r+repmat(U(:,end),[1 size(U0r,2)]);
U0r= fliplr(cat(2,U(:,end),U0r));

w3  = repmat(linspace(0,1,size(U0l,2)),[size(U0l,1) 1]);
w4  = repmat(linspace(1,0,size(U0l,2)),[size(U0l,1) 1]);
U0 = (-w4.*U0l+w3.*(U0r));

V0 = V;