for u = 2:ni
tic
run('swap_matrix_x_720.m');

%rotor culasse meshing 
% First line
r11 = 1;
s11 = 1;
nnum(r11,s11)=(r11-1)*mr+s11;
M(nnum(r11,s11),nnum(r11,s11))      = MP(nnum(r11,s11),nnum(r11,s11))+MP(nnum(r11,s11),nnum(r11,s11)+m-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m);
M(nnum(r11,s11),nnum(r11,s11)+m-1)  = -MP(nnum(r11,s11),nnum(r11,s11)+m-1);
M(nnum(r11,s11),nnum(r11,s11)+1)    = -MP(nnum(r11,s11),nnum(r11,s11)+1);
M(nnum(r11,s11),nnum(r11,s11)+m)    = -MP(nnum(r11,s11),nnum(r11,s11)+m);
for s11 = 2:mr-1
    nnum(r11,s11)=(r11-1)*mr+s11;
    M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11))+MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m);
    M(nnum(r11,s11),nnum(r11,s11)-1)  = -MP(nnum(r11,s11),nnum(r11,s11)-1);
    M(nnum(r11,s11),nnum(r11,s11)+1)  = -MP(nnum(r11,s11),nnum(r11,s11)+1);
    M(nnum(r11,s11),nnum(r11,s11)+m)  = -MP(nnum(r11,s11),nnum(r11,s11)+m);
end
s11 = mr;
nnum(r11,s11)=(r11-1)*mr+s11;
M(nnum(r11,s11),nnum(r11,s11))      = MP(nnum(r11,s11),nnum(r11,s11))+MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m+1);
M(nnum(r11,s11),nnum(r11,s11)-1)    =  -MP(nnum(r11,s11),nnum(r11,s11)-1);
M(nnum(r11,s11),nnum(r11,s11)+m)    =  -MP(nnum(r11,s11),nnum(r11,s11)+m);
M(nnum(r11,s11),nnum(r11,s11)-m+1)  =  -MP(nnum(r11,s11),nnum(r11,s11)-m+1);

% first coulumn of rotor culasse
s11 = 1;
for r11 = 2:pcr-1
    nnum(r11,s11)=(r11-1)*mr+s11;
    M(nnum(r11,s11),nnum(r11,s11))        =  MP(nnum(r11,s11),nnum(r11,s11)+m-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
    M(nnum(r11,s11),nnum(r11,s11)+m-1)    =  -MP(nnum(r11,s11),nnum(r11,s11)+m-1);
    M(nnum(r11,s11),nnum(r11,s11)+1)      =  -MP(nnum(r11,s11),nnum(r11,s11)+1);
    M(nnum(r11,s11),nnum(r11,s11)+m)      =  -MP(nnum(r11,s11),nnum(r11,s11)+m);
    M(nnum(r11,s11),nnum(r11,s11)-m)      =  -MP(nnum(r11,s11),nnum(r11,s11)-m);
end

% last coulumn of rotor culasse
s11 = mr;
for r11 = 2:pcr-1
    nnum(r11,s11)=(r11-1)*mr+s11;
    M(nnum(r11,s11),nnum(r11,s11))        = MP(nnum(r11,s11),nnum(r11,s11)-m+1)+MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
    M(nnum(r11,s11),nnum(r11,s11)-m+1)    =  -MP(nnum(r11,s11),nnum(r11,s11)-m+1);
    M(nnum(r11,s11),nnum(r11,s11)-1)      =  -MP(nnum(r11,s11),nnum(r11,s11)-1);
    M(nnum(r11,s11),nnum(r11,s11)+m)      =  -MP(nnum(r11,s11),nnum(r11,s11)+m);
    M(nnum(r11,s11),nnum(r11,s11)-m)      =  -MP(nnum(r11,s11),nnum(r11,s11)-m);
end

% Rest of elements of rotor culasse
for r11 = 2:pcr-1
    for s11 = 2:mr-1
        nnum(r11,s11)=(r11-1)*mr+s11;
        M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
        M(nnum(r11,s11),nnum(r11,s11)-1) = -MP(nnum(r11,s11),nnum(r11,s11)-1);
        M(nnum(r11,s11),nnum(r11,s11)+1) = -MP(nnum(r11,s11),nnum(r11,s11)+1);
        M(nnum(r11,s11),nnum(r11,s11)+m) = -MP(nnum(r11,s11),nnum(r11,s11)+m);
        M(nnum(r11,s11),nnum(r11,s11)-m) = -MP(nnum(r11,s11),nnum(r11,s11)-m);
    end
end

% Border culasse rotor/teeth rotor
% row number pcr
r11 = pcr;
s11 = 1;
nnum(r11,s11)=(r11-1)*mr+s11;
M(nnum(r11,s11),nnum(r11,s11))      = MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m-1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
M(nnum(r11,s11),nnum(r11,s11)-m)    =  -MP(nnum(r11,s11),nnum(r11,s11)-m);
M(nnum(r11,s11),nnum(r11,s11)+m-1)  =  -MP(nnum(r11,s11),nnum(r11,s11)+m-1);
M(nnum(r11,s11),nnum(r11,s11)+1)    =  -MP(nnum(r11,s11),nnum(r11,s11)+1);
M(nnum(r11,s11),nnum(r11,s11)+m)    =  -MP(nnum(r11,s11),nnum(r11,s11)+m);
for s11 = 2:m0m
    nnum(r11,s11)=(r11-1)*mr+s11;
    M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
    M(nnum(r11,s11),nnum(r11,s11)-1) = -MP(nnum(r11,s11),nnum(r11,s11)-1);
    M(nnum(r11,s11),nnum(r11,s11)+1) = -MP(nnum(r11,s11),nnum(r11,s11)+1);
    M(nnum(r11,s11),nnum(r11,s11)+m) = -MP(nnum(r11,s11),nnum(r11,s11)+m);
    M(nnum(r11,s11),nnum(r11,s11)-m) = -MP(nnum(r11,s11),nnum(r11,s11)-m);
end
for dent = 1:4
    for s11 = (mer+m0m)+1+(dent-1)*(mer+mdr):(mer+m0m+mdr)+(dent-1)*(mer+mdr)
        nnum(r11,s11)=(r11-1)*m+s11;
        M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
        M(nnum(r11,s11),nnum(r11,s11)-1) = -MP(nnum(r11,s11),nnum(r11,s11)-1);
        M(nnum(r11,s11),nnum(r11,s11)+1) = -MP(nnum(r11,s11),nnum(r11,s11)+1);
        M(nnum(r11,s11),nnum(r11,s11)+m) = -MP(nnum(r11,s11),nnum(r11,s11)+m);
        M(nnum(r11,s11),nnum(r11,s11)-m) = -MP(nnum(r11,s11),nnum(r11,s11)-m);
    end
end
for enc = 1:5
    for s11 = m0m+1+(enc-1)*(mer+mdr):(mer+m0m)+(enc-1)*(mer+mdr)
        nnum(r11,s11)=(r11-1)*m+s11;
        M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
        M(nnum(r11,s11),nnum(r11,s11)-1) = -MP(nnum(r11,s11),nnum(r11,s11)-1);
        M(nnum(r11,s11),nnum(r11,s11)+1) = -MP(nnum(r11,s11),nnum(r11,s11)+1);
        M(nnum(r11,s11),nnum(r11,s11)+m) = -MP(nnum(r11,s11),nnum(r11,s11)+m);
        M(nnum(r11,s11),nnum(r11,s11)-m) = -MP(nnum(r11,s11),nnum(r11,s11)-m);
    end
end
for s11 = (5*mer+4*mdr+m0m)+1:mr-1 % here we reached the last half stator tooth
    nnum(r11,s11)=(r11-1)*mr+s11;
    M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
    M(nnum(r11,s11),nnum(r11,s11)-1) = -MP(nnum(r11,s11),nnum(r11,s11)-1);
    M(nnum(r11,s11),nnum(r11,s11)+1) = -MP(nnum(r11,s11),nnum(r11,s11)+1);
    M(nnum(r11,s11),nnum(r11,s11)+m) = -MP(nnum(r11,s11),nnum(r11,s11)+m);
    M(nnum(r11,s11),nnum(r11,s11)-m) = -MP(nnum(r11,s11),nnum(r11,s11)-m);
end
s11 = mr;
nnum(r11,s11)=(r11-1)*mr+s11;
M(nnum(r11,s11),nnum(r11,s11))      = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)-m+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
M(nnum(r11,s11),nnum(r11,s11)-1)   = -MP(nnum(r11,s11),nnum(r11,s11)-1);
M(nnum(r11,s11),nnum(r11,s11)-m+1) = -MP(nnum(r11,s11),nnum(r11,s11)-m+1);
M(nnum(r11,s11),nnum(r11,s11)+m)   = -MP(nnum(r11,s11),nnum(r11,s11)+m);
M(nnum(r11,s11),nnum(r11,s11)-m)   = -MP(nnum(r11,s11),nnum(r11,s11)-m);

% rotor teeth and slots meshing r11 >= pcr+1
% first column s11 = 1
s11 = 1;
for r11 = pcr+1:pr-1
    nnum(r11,s11)=(r11-1)*m+s11;
    M(nnum(r11,s11),nnum(r11,s11))      =  MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m-1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
    M(nnum(r11,s11),nnum(r11,s11)+m-1)  = -MP(nnum(r11,s11),nnum(r11,s11)+m-1);
    M(nnum(r11,s11),nnum(r11,s11)+1)    = -MP(nnum(r11,s11),nnum(r11,s11)+1);
    M(nnum(r11,s11),nnum(r11,s11)+m)    = -MP(nnum(r11,s11),nnum(r11,s11)+m);
    M(nnum(r11,s11),nnum(r11,s11)-m)    = -MP(nnum(r11,s11),nnum(r11,s11)-m);
end
% last row of rotor, airgap border
r11 = pr;
nnum(r11,s11)=(r11-1)*mr+s11;
M(nnum(r11,s11),nnum(r11,s11))      =  MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m-1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
M(nnum(r11,s11),nnum(r11,s11)+m-1)  = -MP(nnum(r11,s11),nnum(r11,s11)+m-1);
M(nnum(r11,s11),nnum(r11,s11)+1)    = -MP(nnum(r11,s11),nnum(r11,s11)+1);
M(nnum(r11,s11),nnum(r11,s11)+m)    = -MP(nnum(r11,s11),nnum(r11,s11)+m);
M(nnum(r11,s11),nnum(r11,s11)-m)    = -MP(nnum(r11,s11),nnum(r11,s11)-m);

r11 = pr+1; % first row of airgap
nnum(r11,s11)=(r11-1)*mr+s11;
M(nnum(r11,s11),nnum(r11,s11))      =  MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m-1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
M(nnum(r11,s11),nnum(r11,s11)+m-1)  = -MP(nnum(r11,s11),nnum(r11,s11)+m-1);
M(nnum(r11,s11),nnum(r11,s11)+1)    = -MP(nnum(r11,s11),nnum(r11,s11)+1);
M(nnum(r11,s11),nnum(r11,s11)+m)    = -MP(nnum(r11,s11),nnum(r11,s11)+m);
M(nnum(r11,s11),nnum(r11,s11)-m)    = -MP(nnum(r11,s11),nnum(r11,s11)-m);

% Last column of rotor (s11 = mr), r11 >= pcr+1
s11 = mr;
for r11 = pcr+1:pr-1
    nnum(r11,s11)=(r11-1)*mr+s11;
    M(nnum(r11,s11),nnum(r11,s11))      = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)-m+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
    M(nnum(r11,s11),nnum(r11,s11)-m+1)  = -MP(nnum(r11,s11),nnum(r11,s11)-m+1);
    M(nnum(r11,s11),nnum(r11,s11)-1)    = -MP(nnum(r11,s11),nnum(r11,s11)-1);
    M(nnum(r11,s11),nnum(r11,s11)+m)    = -MP(nnum(r11,s11),nnum(r11,s11)+m);
    M(nnum(r11,s11),nnum(r11,s11)-m)    = -MP(nnum(r11,s11),nnum(r11,s11)-m);
end
r11 = pr;
nnum(r11,s11)=(r11-1)*mr+s11;
M(nnum(r11,s11),nnum(r11,s11))      = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)-m+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
M(nnum(r11,s11),nnum(r11,s11)-m+1)  = -MP(nnum(r11,s11),nnum(r11,s11)-m+1);
M(nnum(r11,s11),nnum(r11,s11)-1)    = -MP(nnum(r11,s11),nnum(r11,s11)-1);
M(nnum(r11,s11),nnum(r11,s11)+m)    = -MP(nnum(r11,s11),nnum(r11,s11)+m);
M(nnum(r11,s11),nnum(r11,s11)-m)    = -MP(nnum(r11,s11),nnum(r11,s11)-m);

r11 = pr+1;
nnum(r11,s11)=(r11-1)*mr+s11;
M(nnum(r11,s11),nnum(r11,s11))      = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)-m+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
M(nnum(r11,s11),nnum(r11,s11)-m+1)  = -MP(nnum(r11,s11),nnum(r11,s11)-m+1);
M(nnum(r11,s11),nnum(r11,s11)-1)    = -MP(nnum(r11,s11),nnum(r11,s11)-1);
M(nnum(r11,s11),nnum(r11,s11)+m)    = -MP(nnum(r11,s11),nnum(r11,s11)+m);
M(nnum(r11,s11),nnum(r11,s11)-m)    = -MP(nnum(r11,s11),nnum(r11,s11)-m);

% Columns right side borders of rotor teeth r11 >= pcr+1
for dent = 1:5
    s11 = m0m+(dent-1)*(mer+mdr);
    for r11 = pcr+1:pr-1
        nnum(r11,s11)=(r11-1)*m+s11;
        M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
        M(nnum(r11,s11),nnum(r11,s11)-1)    = -MP(nnum(r11,s11),nnum(r11,s11)-1);
        M(nnum(r11,s11),nnum(r11,s11)+1)    = -MP(nnum(r11,s11),nnum(r11,s11)+1);
        M(nnum(r11,s11),nnum(r11,s11)+m)    = -MP(nnum(r11,s11),nnum(r11,s11)+m);
        M(nnum(r11,s11),nnum(r11,s11)-m)    = -MP(nnum(r11,s11),nnum(r11,s11)-m);
    end
    r11 = pr;
    nnum(r11,s11)=(r11-1)*mr+s11;
    M(nnum(r11,s11),nnum(r11,s11))      = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
    M(nnum(r11,s11),nnum(r11,s11)+1)    = -MP(nnum(r11,s11),nnum(r11,s11)+1);
    M(nnum(r11,s11),nnum(r11,s11)-1)    = -MP(nnum(r11,s11),nnum(r11,s11)-1);
    M(nnum(r11,s11),nnum(r11,s11)+m)    = -MP(nnum(r11,s11),nnum(r11,s11)+m);
    M(nnum(r11,s11),nnum(r11,s11)-m)    = -MP(nnum(r11,s11),nnum(r11,s11)-m);
    r11 = pr+1;
    nnum(r11,s11)=(r11-1)*mr+s11;
    M(nnum(r11,s11),nnum(r11,s11))      = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
    M(nnum(r11,s11),nnum(r11,s11)+1)    = -MP(nnum(r11,s11),nnum(r11,s11)+1);
    M(nnum(r11,s11),nnum(r11,s11)-1)    = -MP(nnum(r11,s11),nnum(r11,s11)-1);
    M(nnum(r11,s11),nnum(r11,s11)+m)    = -MP(nnum(r11,s11),nnum(r11,s11)+m);
    M(nnum(r11,s11),nnum(r11,s11)-m)    = -MP(nnum(r11,s11),nnum(r11,s11)-m);
end

% Columns left side borders of teeth  r11 >= pcr+1
for dent = 1:5
    s11 = m0m+mer+1+(dent-1)*(mer+mdr);
    for r11 = pcr+1:pr-1
        nnum(r11,s11)=(r11-1)*m+s11;
        M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
        M(nnum(r11,s11),nnum(r11,s11)-1)  = -MP(nnum(r11,s11),nnum(r11,s11)-1);
        M(nnum(r11,s11),nnum(r11,s11)+1)  = -MP(nnum(r11,s11),nnum(r11,s11)+1);
        M(nnum(r11,s11),nnum(r11,s11)+m)  = -MP(nnum(r11,s11),nnum(r11,s11)+m);
        M(nnum(r11,s11),nnum(r11,s11)-m)  = -MP(nnum(r11,s11),nnum(r11,s11)-m);
    end
    r11 = pr;
    nnum(r11,s11)=(r11-1)*mr+s11;
    M(nnum(r11,s11),nnum(r11,s11))      = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
    M(nnum(r11,s11),nnum(r11,s11)+1)    = -MP(nnum(r11,s11),nnum(r11,s11)+1);
    M(nnum(r11,s11),nnum(r11,s11)-1)    = -MP(nnum(r11,s11),nnum(r11,s11)-1);
    M(nnum(r11,s11),nnum(r11,s11)+m)    = -MP(nnum(r11,s11),nnum(r11,s11)+m);
    M(nnum(r11,s11),nnum(r11,s11)-m)    = -MP(nnum(r11,s11),nnum(r11,s11)-m);
    r11 = pr+1;
    nnum(r11,s11)=(r11-1)*mr+s11;
    M(nnum(r11,s11),nnum(r11,s11))      = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
    M(nnum(r11,s11),nnum(r11,s11)+1)    = -MP(nnum(r11,s11),nnum(r11,s11)+1);
    M(nnum(r11,s11),nnum(r11,s11)-1)    = -MP(nnum(r11,s11),nnum(r11,s11)-1);
    M(nnum(r11,s11),nnum(r11,s11)+m)    = -MP(nnum(r11,s11),nnum(r11,s11)+m);
    M(nnum(r11,s11),nnum(r11,s11)-m)    = -MP(nnum(r11,s11),nnum(r11,s11)-m);
end

% the Half teeth insides without boders
% needed in case m0m > 2
for dent = 1:2
    for s11 = 2+(dent-1)*(5*mer+4*mdr+m0m):m0m-1+(dent-1)*(5*mer+4*mdr+m0m)
        for r11=pcr+1:pr-1
            nnum(r11,s11)=(r11-1)*m+s11;
            M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
            M(nnum(r11,s11),nnum(r11,s11)-1)  = -MP(nnum(r11,s11),nnum(r11,s11)-1);
            M(nnum(r11,s11),nnum(r11,s11)+1)  = -MP(nnum(r11,s11),nnum(r11,s11)+1);
            M(nnum(r11,s11),nnum(r11,s11)+m)  = -MP(nnum(r11,s11),nnum(r11,s11)+m);
            M(nnum(r11,s11),nnum(r11,s11)-m)  = -MP(nnum(r11,s11),nnum(r11,s11)-m);
        end
    r11 = pr;
    nnum(r11,s11)=(r11-1)*m+s11;
    M(nnum(r11,s11),nnum(r11,s11))      = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
    M(nnum(r11,s11),nnum(r11,s11)+1)    = -MP(nnum(r11,s11),nnum(r11,s11)+1);
    M(nnum(r11,s11),nnum(r11,s11)-1)    = -MP(nnum(r11,s11),nnum(r11,s11)-1);
    M(nnum(r11,s11),nnum(r11,s11)+m)    = -MP(nnum(r11,s11),nnum(r11,s11)+m);
    M(nnum(r11,s11),nnum(r11,s11)-m)    = -MP(nnum(r11,s11),nnum(r11,s11)-m);
    r11 = pr+1;
    nnum(r11,s11)=(r11-1)*m+s11;
    M(nnum(r11,s11),nnum(r11,s11))      = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
    M(nnum(r11,s11),nnum(r11,s11)+1)    = -MP(nnum(r11,s11),nnum(r11,s11)+1);
    M(nnum(r11,s11),nnum(r11,s11)-1)    = -MP(nnum(r11,s11),nnum(r11,s11)-1);
    M(nnum(r11,s11),nnum(r11,s11)+m)    = -MP(nnum(r11,s11),nnum(r11,s11)+m);
    M(nnum(r11,s11),nnum(r11,s11)-m)    = -MP(nnum(r11,s11),nnum(r11,s11)-m);
    end
end

% the inside nodes of the teeth
for dent = 1:4
    for s11 = m0m+mer+2+(dent-1)*(mer+mdr):m0m+mer+mdr-1+(dent-1)*(mer+mdr)
        for r11=pcr+1:pr-1
            nnum(r11,s11)=(r11-1)*m+s11;
            M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
            M(nnum(r11,s11),nnum(r11,s11)-1)  = -MP(nnum(r11,s11),nnum(r11,s11)-1);
            M(nnum(r11,s11),nnum(r11,s11)+1)  = -MP(nnum(r11,s11),nnum(r11,s11)+1);
            M(nnum(r11,s11),nnum(r11,s11)+m)  = -MP(nnum(r11,s11),nnum(r11,s11)+m);
            M(nnum(r11,s11),nnum(r11,s11)-m)  = -MP(nnum(r11,s11),nnum(r11,s11)-m);
        end
     r11 = pr;
     nnum(r11,s11)=(r11-1)*m+s11;
     M(nnum(r11,s11),nnum(r11,s11))      = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
     M(nnum(r11,s11),nnum(r11,s11)+1)    = -MP(nnum(r11,s11),nnum(r11,s11)+1);
     M(nnum(r11,s11),nnum(r11,s11)-1)    = -MP(nnum(r11,s11),nnum(r11,s11)-1);
     M(nnum(r11,s11),nnum(r11,s11)+m)    = -MP(nnum(r11,s11),nnum(r11,s11)+m);
     M(nnum(r11,s11),nnum(r11,s11)-m)    = -MP(nnum(r11,s11),nnum(r11,s11)-m);
     r11 = pr+1;
     nnum(r11,s11)=(r11-1)*m+s11;
     M(nnum(r11,s11),nnum(r11,s11))      = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
     M(nnum(r11,s11),nnum(r11,s11)+1)    = -MP(nnum(r11,s11),nnum(r11,s11)+1);
     M(nnum(r11,s11),nnum(r11,s11)-1)    = -MP(nnum(r11,s11),nnum(r11,s11)-1);
     M(nnum(r11,s11),nnum(r11,s11)+m)    = -MP(nnum(r11,s11),nnum(r11,s11)+m);
     M(nnum(r11,s11),nnum(r11,s11)-m)    = -MP(nnum(r11,s11),nnum(r11,s11)-m);
    end
end

% nodes inside slots in contact with the teeth
% left borders of the slots
for enc = 1:5
    s11 = m0m+1+(enc-1)*(mdr+mer);
    r11 = pcr+1;
    nnum(r11,s11)=(r11-1)*m+s11;
    M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
    M(nnum(r11,s11),nnum(r11,s11)-1)  = -MP(nnum(r11,s11),nnum(r11,s11)-1);
    M(nnum(r11,s11),nnum(r11,s11)+1)  = -MP(nnum(r11,s11),nnum(r11,s11)+1);
    M(nnum(r11,s11),nnum(r11,s11)+m)  = -MP(nnum(r11,s11),nnum(r11,s11)+m);
    M(nnum(r11,s11),nnum(r11,s11)-m)  = -MP(nnum(r11,s11),nnum(r11,s11)-m);
    for r11 = pcr+2:pr-1
        nnum(r11,s11)=(r11-1)*m+s11;
        M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
        M(nnum(r11,s11),nnum(r11,s11)-1)  = -MP(nnum(r11,s11),nnum(r11,s11)-1);
        M(nnum(r11,s11),nnum(r11,s11)+1)  = -MP(nnum(r11,s11),nnum(r11,s11)+1);
        M(nnum(r11,s11),nnum(r11,s11)+m)  = -MP(nnum(r11,s11),nnum(r11,s11)+m);
        M(nnum(r11,s11),nnum(r11,s11)-m)  = -MP(nnum(r11,s11),nnum(r11,s11)-m);
    end
    r11 = pr;
    nnum(r11,s11)=(r11-1)*m+s11;
    M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
    M(nnum(r11,s11),nnum(r11,s11)-1)  = -MP(nnum(r11,s11),nnum(r11,s11)-1);
    M(nnum(r11,s11),nnum(r11,s11)+1)  = -MP(nnum(r11,s11),nnum(r11,s11)+1);
    M(nnum(r11,s11),nnum(r11,s11)+m)  = -MP(nnum(r11,s11),nnum(r11,s11)+m);
    M(nnum(r11,s11),nnum(r11,s11)-m)  = -MP(nnum(r11,s11),nnum(r11,s11)-m);
    r11 = pr+1;
    nnum(r11,s11)=(r11-1)*m+s11;
    M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
    M(nnum(r11,s11),nnum(r11,s11)-1)  = -MP(nnum(r11,s11),nnum(r11,s11)-1);
    M(nnum(r11,s11),nnum(r11,s11)+1)  = -MP(nnum(r11,s11),nnum(r11,s11)+1);
    M(nnum(r11,s11),nnum(r11,s11)+m)  = -MP(nnum(r11,s11),nnum(r11,s11)+m);
    M(nnum(r11,s11),nnum(r11,s11)-m)  = -MP(nnum(r11,s11),nnum(r11,s11)-m);
end

% right borders of the slots
for enc = 1:5
    s11 = m0m+mer+(enc-1)*(mer+mdr);
    r11 = pcr+1;
    nnum(r11,s11)=(r11-1)*m+s11;
    M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
    M(nnum(r11,s11),nnum(r11,s11)-1)  = -MP(nnum(r11,s11),nnum(r11,s11)-1);
    M(nnum(r11,s11),nnum(r11,s11)+1)  = -MP(nnum(r11,s11),nnum(r11,s11)+1);
    M(nnum(r11,s11),nnum(r11,s11)+m)  = -MP(nnum(r11,s11),nnum(r11,s11)+m);
    M(nnum(r11,s11),nnum(r11,s11)-m)  = -MP(nnum(r11,s11),nnum(r11,s11)-m);
    for r11 = pcr+2:pr-1
        nnum(r11,s11)=(r11-1)*m+s11;
        M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
        M(nnum(r11,s11),nnum(r11,s11)-1)  = -MP(nnum(r11,s11),nnum(r11,s11)-1);
        M(nnum(r11,s11),nnum(r11,s11)+1)  = -MP(nnum(r11,s11),nnum(r11,s11)+1);
        M(nnum(r11,s11),nnum(r11,s11)+m)  = -MP(nnum(r11,s11),nnum(r11,s11)+m);
        M(nnum(r11,s11),nnum(r11,s11)-m)  = -MP(nnum(r11,s11),nnum(r11,s11)-m);
    end
    r11 = pr;
    nnum(r11,s11)=(r11-1)*m+s11;
    M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
    M(nnum(r11,s11),nnum(r11,s11)-1)  = -MP(nnum(r11,s11),nnum(r11,s11)-1);
    M(nnum(r11,s11),nnum(r11,s11)+1)  = -MP(nnum(r11,s11),nnum(r11,s11)+1);
    M(nnum(r11,s11),nnum(r11,s11)+m)  = -MP(nnum(r11,s11),nnum(r11,s11)+m);
    M(nnum(r11,s11),nnum(r11,s11)-m)  = -MP(nnum(r11,s11),nnum(r11,s11)-m);
    r11 = pr+1;
    nnum(r11,s11)=(r11-1)*m+s11;
    M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
    M(nnum(r11,s11),nnum(r11,s11)-1)  = -MP(nnum(r11,s11),nnum(r11,s11)-1);
    M(nnum(r11,s11),nnum(r11,s11)+1)  = -MP(nnum(r11,s11),nnum(r11,s11)+1);
    M(nnum(r11,s11),nnum(r11,s11)+m)  = -MP(nnum(r11,s11),nnum(r11,s11)+m);
    M(nnum(r11,s11),nnum(r11,s11)-m)  = -MP(nnum(r11,s11),nnum(r11,s11)-m);
    
end

% inner nodes of the rotor slots
for enc = 1:5
    for s11 = m0m+2+(enc-1)*(mer+mdr):m0m+mer-1+(enc-1)*(mer+mdr)
        r11 = pcr+1;
        nnum(r11,s11)=(r11-1)*m+s11;
        M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
        M(nnum(r11,s11),nnum(r11,s11)-1)  = -MP(nnum(r11,s11),nnum(r11,s11)-1);
        M(nnum(r11,s11),nnum(r11,s11)+1)  = -MP(nnum(r11,s11),nnum(r11,s11)+1);
        M(nnum(r11,s11),nnum(r11,s11)+m)  = -MP(nnum(r11,s11),nnum(r11,s11)+m);
        M(nnum(r11,s11),nnum(r11,s11)-m)  = -MP(nnum(r11,s11),nnum(r11,s11)-m);
        for r11 = pcr+2:pr-1
            nnum(r11,s11)=(r11-1)*m+s11;
            M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
            M(nnum(r11,s11),nnum(r11,s11)-1)  = -MP(nnum(r11,s11),nnum(r11,s11)-1);
            M(nnum(r11,s11),nnum(r11,s11)+1)  = -MP(nnum(r11,s11),nnum(r11,s11)+1);
            M(nnum(r11,s11),nnum(r11,s11)+m)  = -MP(nnum(r11,s11),nnum(r11,s11)+m);
            M(nnum(r11,s11),nnum(r11,s11)-m)  = -MP(nnum(r11,s11),nnum(r11,s11)-m);
        end
        r11 = pr;
        nnum(r11,s11)=(r11-1)*m+s11;
        M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
        M(nnum(r11,s11),nnum(r11,s11)-1)  = -MP(nnum(r11,s11),nnum(r11,s11)-1);
        M(nnum(r11,s11),nnum(r11,s11)+1)  = -MP(nnum(r11,s11),nnum(r11,s11)+1);
        M(nnum(r11,s11),nnum(r11,s11)+m)  = -MP(nnum(r11,s11),nnum(r11,s11)+m);
        M(nnum(r11,s11),nnum(r11,s11)-m)  = -MP(nnum(r11,s11),nnum(r11,s11)-m);
        r11 = pr+1;
        nnum(r11,s11)=(r11-1)*m+s11;
        M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
        M(nnum(r11,s11),nnum(r11,s11)-1)  = -MP(nnum(r11,s11),nnum(r11,s11)-1);
        M(nnum(r11,s11),nnum(r11,s11)+1)  = -MP(nnum(r11,s11),nnum(r11,s11)+1);
        M(nnum(r11,s11),nnum(r11,s11)+m)  = -MP(nnum(r11,s11),nnum(r11,s11)+m);
        M(nnum(r11,s11),nnum(r11,s11)-m)  = -MP(nnum(r11,s11),nnum(r11,s11)-m);
    end
end


C=M\E';

Ups = zeros(p+p1,m);
nnum3 = zeros(p+p1,m);

for i=1:(p+p1)
    for j=1:m
        nnum3(i,j)=(i-1)*m+j;
        Ups(i,j) = C(nnum3(i,j));        
    end
end
Phi = zeros(length(rangeTooth1),12);
i=1;

% tooth 1
for s11 = rangeTooth1
    nnum11=(p1+5)*m+s11;
    nnum12=(p1+5+1)*m+s11;
    Phi(i,1)=(Ups((p1+5),s11)-Ups((p1+5+1),s11))*MP(nnum11,nnum12);
    i = i+1;
end
i=1;
% tooth 2
for s11 = rangeTooth2
    nnum21=(p1+5)*m+s11;
    nnum22=(p1+5+1)*m+s11;
    Phi(i,2)=(Ups((p1+5),s11)-Ups((p1+5+1),s11))*MP(nnum21,nnum22);
    i=i+1;
end
i=1;
% tooth 3
for s11 = rangeTooth3
    nnum31=(p1+5)*m+s11;
    nnum32=(p1+5+1)*m+s11;
    Phi(i,3)=(Ups((p1+5),s11)-Ups((p1+5+1),s11))*MP(nnum31,nnum32);
    i=i+1;
end
i=1;
% tooth 4
for s11 = rangeTooth4
    nnum41=(p1+5)*m+s11;
    nnum42=(p1+5+1)*m+s11;
    Phi(i,4)=(Ups((p1+5),s11)-Ups((p1+5+1),s11))*MP(nnum41,nnum42);
    i=i+1;
end
i=1;
%tooth 5
for s11 = rangeTooth5
    nnum51=(p1+5)*m+s11;
    nnum52=(p1+5+1)*m+s11;
    Phi(i,5)=(Ups((p1+5),s11)-Ups((p1+5+1),s11))*MP(nnum51,nnum52);
    i=i+1;
end
i=1;
% tooth 6
for s11 = rangeTooth6
    nnum61=(p1+5)*m+s11;
    nnum62=(p1+5+1)*m+s11;
    Phi(i,6)=(Ups((p1+5),s11)-Ups((p1+5+1),s11))*MP(nnum61,nnum62);
    i=i+1;
end
i=1;
% tooth 7
for s11 = rangeTooth7
    nnum71=(p1+5)*m+s11;
    nnum72=(p1+5+1)*m+s11;
    Phi(i,7)=(Ups((p1+5),s11)-Ups((p1+5+1),s11))*MP(nnum71,nnum72);
    i=i+1;
end
i=1;
% tooth 8
for s11 = rangeTooth8
    nnum81=(p1+5)*m+s11;
    nnum82=(p1+5+1)*m+s11;
    Phi(i,8)=(Ups((p1+5),s11)-Ups((p1+5+1),s11))*MP(nnum81,nnum82);
    i=i+1;
end
i=1;
% tooth 9
for s11 = rangeTooth9
    nnum91=(p1+5)*m+s11;
    nnum92=(p1+5+1)*m+s11;
    Phi(i,9)=(Ups((p1+5),s11)-Ups((p1+5+1),s11))*MP(nnum91,nnum92);
    i=i+1;
end
i=1;
% tooth 10
for s11 = rangeTooth10
    nnum101=(p1+5)*m+s11;
    nnum102=(p1+5+1)*m+s11;
    Phi(i,10)=(Ups((p1+5),s11)-Ups((p1+5+1),s11))*MP(nnum101,nnum102);
    i=i+1;
end
i=1;
% tooth 11
for s11 = rangeTooth11
    nnum1111=(p1+5)*m+s11;
    nnum1112=(p1+5+1)*m+s11;
    Phi(i,11)=(Ups((p1+5),s11)-Ups((p1+5+1),s11))*MP(nnum1111,nnum1112);
    i=i+1;
end
i=1;
% tooth 12
for s11 = rangeTooth12
    nnum1121=(p1+5)*m+s11;
    nnum1122=(p1+5+1)*m+s11;
    Phi(i,12)=(Ups((p1+5),s11)-Ups((p1+5+1),s11))*MP(nnum1121,nnum1122);
    i=i+1;
end


% Calculation of magnetic flux density in the air gap
for s11 = 2:m-1
    nnumB1   = (p1-1)*m+s11;
    Bx(s11,u) = mu0*(C(nnumB1-1)-C(nnumB1+1))/2/l;
end

nnumB1   = (p1-1)*m+1;
Bx(1,u) = mu0*(-C(nnumB1-1+m)-C(nnumB1+1))/2/l;
Bx(m,u) = mu0*(C(nnumB1-1+m-1)+C(nnumB1))/2/l;


for i = 1:m
    nnumB31   = (p1)*m+i;
    nnumB22   = (p1-1)*m+i;
    nnumB23   = (p1-2)*m+i;
    nnumB24   = (p1-3)*m+i;
    By1(i,u) = (C(nnumB23)-C(nnumB22))*MP(nnumB23,nnumB22)/(l*la);
    By2(i,u) = (C(nnumB22)-C(nnumB31))*MP(nnumB22,nnumB31)/(l*la);
    By3(i,u) = (C(nnumB24)-C(nnumB23))*MP(nnumB24,nnumB23)/(l*la);
    By(i,u)  = (By1(i,u)+By3(i,u)+By2(i,u))/3;
end

FluxB(u)     = sum(Phi(:,2))+sum(Phi(:,3))+sum(Phi(:,8))+sum(Phi(:,9));
FluxA(u)     = sum(Phi(:,1))+sum(Phi(:,12))+sum(Phi(:,6))+sum(Phi(:,7));
FluxC(u)     = sum(Phi(:,4))+sum(Phi(:,5))+sum(Phi(:,10))+sum(Phi(:,11));
toc
fprintf('iteration number: %d\n',u);
end 
figure
plot(1:u,FluxA,1:u,FluxB,1:u,FluxC);