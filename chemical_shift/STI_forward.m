function STI_forward(filename)
    % Compute field shift delta
    % Input:
    %       chi (3x3 'cell'), is the phantom generated by `gen_phantom.m`
    % Save:
    %           % STIParams---
    %           OriNum: Orientation Number
    %           fov, sizeVol, VoxSize: field of view, size of Volume, voxel size
    %           filename_deltaB: cell array containing the filenames of deltaB maps
    %           maskMSA: White Matter Mask
    %           QSM:
    mat = load(['data/',filename]);
    Params = mat.Params;
    chi = mat.chi;
    SNR = 30;
    noiselevel = 0.033;
    Params.noiselevel = noiselevel;
    Params.Rorder = [1, 2, 3];     % rotation order   
    Params.shift = -3.5;           % chemical shift -3.5 ppm
    BrainMask = mat.BrainMask;

    OriAngle = '20-40';
    OriAngleRLarray = [0, -20, -20, -20, 40, 40, 40];
    OriAngleAParray = [0,   0,   0,   0,  0,  0,  0];
    OriAngleFHarray = [0,   0, 120, 240,  0,120,240];
    STIParams.filename_deltaB = cell(length(OriAngleRLarray),1);
    for OrientInd = 1:length(OriAngleRLarray)

        Params.AngRL = OriAngleRLarray(OrientInd);
        Params.AngAP = OriAngleAParray(OrientInd);
        Params.AngFH = OriAngleFHarray(OrientInd);
        Params.TAng = Rmatrix_arb(Params.AngRL, Params.AngAP, Params.AngFH, Params.Rorder);     % rotation matrix

        STIParams.filename_deltaB{OrientInd} = fullfile(Params.pathname, ['O', num2str(OrientInd), '_SNR', num2str(SNR), '_Ang', OriAngle, '.mat']);
        Params.filename_deltaB = STIParams.filename_deltaB{OrientInd};
        
        % simulate field map
        deltaB = chitensor2delta(chi, Params);
        % add chemical shift
        deltaB = deltaB + (chi{1,1}~=0).*Param.shift;
        % add noise
        deltaB = deltaB + randn(size(chi{1,1}))*Params.noiselevel.*(std(deltaB(:)));        
        
        save(Params.filename_deltaB, 'BrainMask', 'deltaB', 'Params')
        disp([Params.filename_deltaB, ' saved.'])

    end

    STIParams.OriNum = length(OriAngleRLarray);
    STIParams.fov = Params.fov;
    STIParams.sizeVol = Params.sizeVol;
    STIParams.VoxSize = Params.voxSize;
    STIParams.maskMSA = mat.chiani~=0;
    STIParams.QSM = mat.chiavg;
    save('data/STIParams.mat', 'STIParams');
    disp('data/STIParams.mat saved');
end