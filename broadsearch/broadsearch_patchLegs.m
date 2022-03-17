function dataTrajUF = broadsearch_patchLegs(dataTrajUF, leg1, leg2, patchName,Jmax)
%BROADSEARCH_PATCHLEGS Patch two leg sets of Lambert solutions
%
%   Assumptions/Warnings:
%   	[NOTE: This function is not meant to be standalone]
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Dependencies:
%       1. broadsearch_patchLegs()
%       2. getPlanetName()
%       3. bplanefromVi1Vi2()
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%

    global pcd

    leg1 = dataTrajUF.(leg1);
    leg2 = dataTrajUF.(leg2);
    pnum = str2double(patchName(2:end));

    dataTrajUF.([patchName]).db = leg1.db;
    dataTrajUF.([patchName]).ab = leg2.ab;
    dataTrajUF.([patchName]).fb = leg1.ab;

    if pnum==1
        cost(1:height(leg1.legdata),1)     = NaN;
        refleg(1:height(leg1.legdata),1:2) = NaN;
        flyby(1:height(leg1.legdata),1:4)  = NaN;
    elseif pnum>1

        cost   = dataTrajUF.(['p',num2str(pnum-1)]).cost;
        cost(:,((pnum-1)*1)+1) = NaN;

        refleg = dataTrajUF.(['p',num2str(pnum-1)]).refleg;
        refleg(:,((pnum-1)*2)+1:((pnum-1)*2)+2) = NaN;

        flyby  = dataTrajUF.(['p',num2str(pnum-1)]).flyby;
        flyby(:,((pnum-1)*4)+1:((pnum-1)*4)+4) = NaN;

    end

    for i=1:height(leg1.legdata)

        for j=1:height(leg2.legdata)

            % Match Dates
            date1 = leg1.legdata(i,3);
            date2 = leg2.legdata(j,2);

            if date1==date2
                % Get Cost
                J1 = leg1.legdata(i,6);
                J2 = leg2.legdata(j,5);
                J  = abs(J2-J1);


                if J<=Jmax

                    % Compute Flyby Parameters
                    fbb  = leg1.ab;
                    vi1_ = leg1.legdata(i,10:12).';
                    vi2_ = leg2.legdata(j,7:9).';   
                    mufb = pcd.(getPlanetName(fbb)).mu;
                    bpdata = bplanefromVi1Vi2(mufb, vi1_, vi2_);

                    % Close Approach Radius and Turning Angle
                    rp = bpdata.rp;
                    d  = bpdata.d;

                    if rp > pcd.(getPlanetName(fbb)).r
                        legdata(i,1)     = i;
                        legdata(i,2)     = leg1.legdata(i,2);
                        legdata(i,3)     = leg2.legdata(j,3);
                        legdata(i,4)     = (legdata(i,3) - legdata(i,2))*86400;
                        legdata(i,5)     = leg1.legdata(i,5);
                        legdata(i,6)     = leg2.legdata(j,6);
                        legdata(i,7:9)   = leg1.legdata(i,7:9);
                        legdata(i,10:12) = leg2.legdata(j,10:12);
                        cost(i,((pnum-1)*1)+1)  = J;
                        refleg(i,((pnum-1)*2)+1) = i;
                        refleg(i,((pnum-1)*2)+1+1) = j;
                        flyby(i,((pnum-1)*4)+1) = date1;
                        flyby(i,((pnum-1)*4)+2) = rp;
                        flyby(i,((pnum-1)*4)+3) = rp-pcd.(getPlanetName(fbb)).r;
                        flyby(i,((pnum-1)*4)+4) = d;
                    else
                        legdata(i,1:12) = NaN;
                        cost(i,((pnum-1)*1)+1)  = NaN;
                        refleg(i,((pnum-1)*2)+1) = NaN;
                        refleg(i,((pnum-1)*2)+1+1) = NaN;
                        flyby(i,((pnum-1)*3)+1) = NaN;
                        flyby(i,((pnum-1)*3)+1+1) = NaN;
                        flyby(i,((pnum-1)*3)+1+2) = NaN;
                    end
                end
            end
        end
    end

    if isnan(cost(:,((pnum-1)*1)+1))
        dataTrajUF.([patchName]).error = true;
    else
        dataTrajUF.([patchName]).error   = false;
        dataTrajUF.([patchName]).legdata = legdata;
        dataTrajUF.([patchName]).flyby   = flyby;
        dataTrajUF.([patchName]).refleg  = refleg;
        dataTrajUF.([patchName]).cost    = cost;
    end

end