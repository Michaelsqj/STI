function STI_forward(Params)
    SNR = 30;
    noiselevel = 0.033;
    Params.noiselevel = noiselevel;
    Params.Rorder = [1, 2, 3];     % rotation order    

    OriAngle = '20-40';
    OriAngleRLarray = [0, -20, -20, -20, 40, 40, 40];
    OriAngleAParray = [0,   0,   0,   0,  0,  0,  0];
    OriAngleFHarray = [0,   0, 120, 240,  0,120,240];

    H0Lab = [0, 0, 1]';
    H0Sub = zeros(3, length(OriAngleRLarray));

    for OrientInd = 1:length(OriAngleRLarray)

        Params.AngRL = OriAngleRLarray(OrientInd);
        Params.AngAP = OriAngleAParray(OrientInd);
        Params.AngFH = OriAngleFHarray(OrientInd);
        Params.TAng = Rmatrix_arb(Params.AngRL, Params.AngAP, Params.AngFH, Params.Rorder);     % rotation matrix

        H0Sub(:,OrientInd) = Params.TAng'*H0Lab;
        Params.filename_deltaB = fullfile(Params.pathname, ['O', num2str(OrientInd), '_SNR', num2str(SNR), '_Ang', OriAngle, '.mat']);

        % simulate field map
        deltaB = chitensor2delta(chi11, chi12, chi13, chi22, chi23, chi33, Params);
        deltaB = deltaB + randn(size(chi11))*Params.noiselevel.*(std(deltaB(:)));

        save(Params.filename_deltaB, 'BrainMask', 'deltaB', 'Params')
        disp([Params.filename_deltaB, ' saved.'])

    end
end

