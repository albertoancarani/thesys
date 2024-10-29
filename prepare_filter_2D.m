
function [Hs, H] = prepare_filter_2D(rmin, nelx, nely)
    nele = nelx * nely;
    iH = ones(nele * (2 * (ceil(rmin) - 1) + 1)^2, 1);
    jH = ones(size(iH));
    sH = zeros(size(iH));
    k = 0;
    
    for i1 = 1:nelx
        for j1 = 1:nely
            e1 = (i1 - 1) * nely + j1;
            for i2 = max(i1 - (ceil(rmin) - 1), 1):min(i1 + (ceil(rmin) - 1), nelx)
                for j2 = max(j1 - (ceil(rmin) - 1), 1):min(j1 + (ceil(rmin) - 1), nely)
                    e2 = (i2 - 1) * nely + j2;
                    k = k + 1;
                    iH(k) = e1;
                    jH(k) = e2;
                    dist = sqrt((i1 - i2)^2 + (j1 - j2)^2);
                    sH(k) = max(0, rmin - dist);
                end
            end
        end
    end
    
    
    H = sparse(iH, jH, sH);

    
    Hs = sum(H, 2);
    
    if any(Hs == 0)
        warning('Hs contiene zeri.');
        Hs(Hs == 0) = 1e-6;  
    end
end