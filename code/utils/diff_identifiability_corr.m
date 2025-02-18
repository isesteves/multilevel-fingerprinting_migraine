function corr_data = diff_identifiability_corr(corr_matrix, template)
    
    % template is coded as natural numbers 
    [N, edges] = histcounts(template(:), [unique(template)', max(unique(template))+1]);
    cases = edges(unique(template)>=1);
    cases_length = N(unique(template)>=1);
    nr_cases = length(cases);
    maxlength = max(cases_length);
    
    corr_data = NaN*ones(maxlength,nr_cases);
    for k = 1:nr_cases
        corr_data(1:cases_length(k), k) = corr_matrix(template == k);       
    end
end
    