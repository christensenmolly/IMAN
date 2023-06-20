def inten_along_ellipse(data, mask, params, r):
    xCen = params["xCen"]
    yCen = params["yCen"]
    ell = params["ell"]
    posang = params["pa"]
    ySize, xSize = data.shape
    cospa = cos(posang)
    sinpa = sin(posang)
    q = 1.0 - ell
    rq = r * q
    eList = linspace(0, 2*pi, 2*r)
    inten = zeros_like(eList)
    for idx, E in enumerate(eList):
        cose = cos(E)
        sine = sin(E)
        fx, ix = modf(xCen + r*cose*cospa - rq*sine*sinpa)
        fy, iy = modf(yCen + rq*sine*cospa + r*cose*sinpa)
        iix = int(ix)
        iiy = int(iy)
        if (iix>0) and (iix<xSize-1) and (iiy>0) and (iiy<ySize-1) and (mask[iiy][iix]==0):
            # supbixel averaging
            inten[idx] = ((1.0-fy)*(1.0-fx)*data[iiy][iix] + 
                          fy*(1.0-fx)*data[iiy+1][iix] + 
                          fx*(1.0-fy)*data[iiy][iix+1] + 
                          fy*fx*data[iiy+1][iix+1])
    return inten
