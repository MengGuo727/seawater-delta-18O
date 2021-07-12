function delta_18_all = phi_to_delta(M_18, M, R_smow, type)

nt = length(M_18);

phi_all = nan(1,nt);
R_18_all = nan(1,nt);
delta_18_all = nan(1,nt);

if (type ==1)
    for i = 1:nt
        phi_all(i) = M_18(i) / M(i);
        R_18_all(i) = phi_all(i) / (1 - phi_all(i));
        delta_18_all(i) = (R_18_all(i) / R_smow - 1) * 1000;
    end
else
    for i = 1:nt
        phi_all(i) = M_18(i) / M;
        R_18_all(i) = phi_all(i) / (1 - phi_all(i));
        delta_18_all(i) = (R_18_all(i) / R_smow - 1) * 1000;
    end
end
