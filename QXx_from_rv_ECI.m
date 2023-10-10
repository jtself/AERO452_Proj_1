function QXx = QXx_from_rv_ECI(r,v)

% outputs DCM from input r, v vectors in ECI

% convert from ECI to LVLH
    % CALCULATE ANGULAR MOMENTUM VECTOR OF S/C A
    hA = cross(r,v); % km2/s

    % CALCULATE UNIT VECTORS OF LVLH FRAME
    i_hat = r / norm(r); 
    k_hat = h / norm(h); % mag of h never changes (assume no perturbations)
    j_hat = cross(k_hat,i_hat);
    
    % CALCULATE ORTHOGONAL DIRECTION COSINE MATRIX (QXx) using Eq. 7.11
    QXx = [i_hat'; j_hat'; k_hat'];

end