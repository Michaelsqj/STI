function [RE_MMS, RE_MSA, AE] = evaluate(chi_r, chi_t)
    % Compute relative error and angular error
    % RE (Relative error) = sqrt(sum((chi_t - chi_r).^2)/sum(chi_t.^2))
    %                                 chi can be MMS or MSA
    % AE (Angular error) = arccos(v_t - v_r)   v is principle eig vector
    %
    % Input:
    %       chi_r: reconstructed chi, class cell, 3x3
    %       chi_t: ground truth chi, class cell, 3x3

    [chi1_r,chi2_r, chi3_r, eig1_r, ~, ~] = tensor2eig(chi_r);
    [chi1_t,chi2_t, chi3_t, eig1_t, ~, ~] = tensor2eig(chi_t);
    
    chiavg_r = (chi1_r+chi_2_r+chi3_r)/3;
    chiavg_t = (chi1_t+chi_2_t+chi3_t)/3;
    chiani_r = chi1_r-(chi2_r+chi3_r)/2;
    chiani_t = chi1_t-(chi2_t+chi3_t)/2;
    % ---------------- Compute RE ----------------
    RE_MMS = sqrt(sum((chiavg_r-chiavg_t).^2,'all')/sum(chiavg_t.^2,'all'));
    RE_MSA = sqrt(sum((chiani_r-chiani_t).^2,'all')/sum(chiani_t.^2,'all'));
    
    % ---------------- Compute AE ----------------
    eig1_r = reshape(eig1_r, 3, numel(eig1_r)/3);
    eig1_t = reshape(eig1_t, 3, numel(eig1_t)/3);
    AE = 0;
    for i = 1:numel(eig1_r)/3
        AE = AE + acosd(eig1_r, eig1_t);
    end
    AE = AE/(numel(eig1_r)/3);
end