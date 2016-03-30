function [eval] = laplace_outlier_detect(eval,dx,dy,passes)

[S1,S2] = size(eval);

for n = 1:passes
    L    = socdiff2(eval,dy,1)+socdiff2(eval,dx,1);
    Eval = zeros(size(L));
    
    q = [1 1];
    
    for i = 1:size(L,1)
        for j = 1:size(L,2)
            if Eval(i,j)==0
                Imin = max([i-q(2) 1   ]);
                Imax = min([i+q(2) S1]);
                
                Jmin = max([j-q(1) 1   ]);
                Jmax = min([j+q(1) S2]);
                
                Iind = Imin:Imax;
                
                Jind = Jmin:Jmax;
                
                block = L(Iind,Jind);
                block(q+1,q+1) = nan;
            end
                        
            Ipos = find(Iind==i);
            Jpos = find(Jind==j);
            
            block_mean = nanmean(block(:));
%             block_std  = nanstd( block(:));
            block_std  = sqrt(nanmean(block(:).^2));
            
%             if L(i,j) > block_mean+block_std || L(i,j) < block_mean-block_std
%                 Eval(i,j) = 1;
%             end 
            if L(i,j) > 2*block_std || L(i,j) < -2*block_std
                Eval(i,j) = 1;
            end
        end
    end
    
    evaltemp = eval;
        
    for i = 1:size(L,1)
        for j = 1:size(L,2)
            if Eval(i,j) == 1
                Imin = max([i-q(2) 1   ]);
                Imax = min([i+q(2) S1]);
                
                Jmin = max([j-q(1) 1   ]);
                Jmax = min([j+q(1) S2]);
                
                Iind = Imin:Imax;
                
                Jind = Jmin:Jmax;
                
                block = eval(Iind,Jind);
                block(q+1,q+1) = nan;
                
                blockind = Eval(Iind,Jind);
                
                block(blockind == 1) = nan;
                                
                evaltemp(i,j) = nanmean(block(:));
                
                if isnan(evaltemp(i,j))
                    evaltemp(i,j) = eval(i,j);
                end
            end
        end
    end
    eval = evaltemp;
end