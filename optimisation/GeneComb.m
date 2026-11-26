function [result_comb] = GeneComb(model,FVA_mut,r_up,r_null,up,result_comb,x)


DM_IgG = findRxnIDs(model,'DM_igg[g]');

b(:,1) = model.lb;
b(:,2) = model.ub;


b_temp = b;
b_temp(r_up,1) = b(r_up,1) + (b(r_up,2)-b(r_up,1)) * up;
b_temp(r_null,:) = 0;
    
model.lb = b_temp(:,1);
model.ub = b_temp(:,2);
sol_init = optimizeCbModel(model);
bio_init = sol_init.f;
dm_init = sol_init.x(DM_IgG);


% Number of combinations

total_reg = size(r_up,2) + size(r_null,2);
% result_comb = zeros(total_reg);

found = 0;
js = 2;

while ( js <= total_reg && found == 0 )
    
    comb_c = js;

    t = [r_null r_up];
    combinations = nchoosek(t,comb_c);

    igg_id = findRxnIDs(model,'DM_igg[g]');

    % result_comb = zeros(comb_c);
%     counter = 0;

    
    i = 1;
    size(combinations,1)
    while (i <= size(combinations,1) && found == 0 && i < 3000 )       
        FVA_mut_temp = FVA_mut;
        
        for j = 1 : size(combinations(i,:),2)
            if find(ismember(r_null,combinations(i,j)))
                FVA_mut_temp( combinations(i,j),:) = 0;
            else
                FVA_mut_temp( combinations(i,j),1) = FVA_mut( combinations(i,j),1) + (FVA_mut( combinations(i,j),2)-FVA_mut( combinations(i,j),1)) * up;
            end
        end
        
        model.lb = FVA_mut_temp(:,1);
        model.ub = FVA_mut_temp(:,2);    
        sol = optimizeCbModel(model);
        igg = sol.x(igg_id);
        
        if igg > 0.1 * dm_init
%             counter = counter + 1;
%             result_comb(end + 1,:) = [combinations(i,:), js];
            result_comb(end + 1,:) = [combinations(i,:)];
            found = 1;
        end
        i = i + 1;
    end    
    js = js + 1;
end
result_comb;

end

