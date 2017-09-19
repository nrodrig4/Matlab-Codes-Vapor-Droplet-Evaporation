function [mu,compound,dirext,M,C_v,R] = CompoundInfo(compoundValues)
%DataBase of Compound Information
   
C = compoundValues;

    if C==1
        mu=507093.5123386746; %unit?
        compound='methanol';
        dirext = 'methanol3200_2500';
        M = 32.04; %g/mol
        C_v = 0.2; %unit?
        R = 6.5; %mm
    elseif C==2
        mu=3.80213e5;
        compound='methanol';
        dirext = 'methanol1150_800';
        M = 32.04; %g/mol
    elseif C==3
        mu=1.954e6;
        compound='hexane';
        dirext = 'hexane3100_2800';
    elseif C==4
        mu=1.108e5;
        compound='hexane';
        dirext = 'hexane1520_1420';
    elseif C==5
        mu=1.829e6;
        compound='3MP';
        dirext = '3MP3100_2800';
    elseif C==6
        mu=1.356e5;
        compound='3MP';
        dirext = '3MP1520_1420';
    end




end

