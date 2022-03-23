function out = getGMATmodJD(jd)
%GETGMATMODJD Returns Modified Julian Data Compatible w/ GMAT
    out = jd-2430000.0;
end