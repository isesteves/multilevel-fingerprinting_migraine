function x_norm = normalize_zero_one(x) 

x_norm = (x-min(x(:)))./(max(x(:))-min(x(:)));