function y = add_noise(X, SNR)
    % add gaussian noise per SNR(dB)
    sizeVol = size(X);
    X = X(:);
    L = length(X);
    SNR_lin = 10^(SNR/10);
    Px = sum(X.^2)/L;
    noise = sqrt(Px/SNR_lin).*randn(L,1);
    y = X + noise;
    y = reshape(y, sizeVol(1), sizeVol(2), sizeVol(3));
end