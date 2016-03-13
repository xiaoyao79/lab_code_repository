function [Uval,Vval] = vel_replace(Uval,Vval,Evalval,winsize)

maxSearch = floor( (max(winsize(:))-1)/2 );

S = size(Uval);

[X,Y] = meshgrid(1:S(2),1:S(1));

%replacement
for i=1:S(1)
    for j=1:S(2)
        if Evalval(i,j)>0
            %initialize replacement search size
            q=0;
            s=0;

            %get replacement block with at least 8 valid points, but stop if region grows bigger than UOD test region
            while s==0
                q=q+1;
                Imin = max([i-q 1   ]);
                Imax = min([i+q S(1)]);
                Jmin = max([j-q 1   ]);
                Jmax = min([j+q S(2)]);
                Iind = Imin:Imax;
                Jind = Jmin:Jmax;
                Ublock = Uval(Iind,Jind);
                if q >= maxSearch || length(Ublock(~isnan(Ublock)))>=8 
                    Xblock = X(Iind,Jind)-X(i,j);
                    Yblock = Y(Iind,Jind)-Y(i,j);
                    Vblock = Vval(Iind,Jind);
                    s=1;
                end
            end
            
            %distance from erroneous vector
            Dblock = (Xblock.^2+Yblock.^2).^-0.5;
            Dblock(isnan(Ublock))=nan;
            Dblock(isinf(Dblock))=nan;

            %validated vector
%             Uval(i,j) = nansum(nansum(Dblock.*Ublock))/nansum(nansum(Dblock));
%             Vval(i,j) = nansum(nansum(Dblock.*Vblock))/nansum(nansum(Dblock));
            
            Uval(i,j) = nanmedian(Dblock(:).*Ublock(:));
            Vval(i,j) = nanmedian(Dblock(:).*Vblock(:));
        end
    end
end