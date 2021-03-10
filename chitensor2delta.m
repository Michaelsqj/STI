function delta = chitensor2delta()
    chi11k = fftn(chi11);                     % tensor components in k space
    chi12k = fftn(chi12);
    chi13k = fftn(chi13);
    chi22k = fftn(chi22);
    chi23k = fftn(chi23);
    chi33k = fftn(chi33);

    [a11, a12, a13, a22, a23, a33] = conv_kernel_tensor(Params);

    delta = a11.*chi11k + a12.*chi12k + a13.*chi13k + ...
                a22.*chi22k + a23.*chi23k + a33.*chi33k ;

    delta = real(ifftn(delta));
end