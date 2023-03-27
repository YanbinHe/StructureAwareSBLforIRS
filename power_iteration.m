function [ld,lv,rd] = power_iteration(mat)

[m,~] = size(mat);
ld = ones(m,1);
ld_old = zeros(m,1);
while norm(ld_old - ld) > 1e-4
    ld_old = ld;
    
    ld = mat*(mat'*ld);
    ld = ld/norm(ld);
end
rdlv = mat'*ld;
lv = norm(rdlv);
rd = rdlv/lv;
end

