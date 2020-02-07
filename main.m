tic
clearvars

Rint = 133e-03;
Rext = 186e-03;
R = (Rint+Rext)/2;   % mean radius

ni = 300;            % number of iterations

awds  = 8.75*pi/180; % stator tooth angle in radians 
awes = 6.25*pi/180;  % stator slot angle in radians
awdr = 15*pi/180;    % rotor tooth angle in radians
awer = 21*pi/180;    % rotor slot angle in radians

wds = awds*R;     % stator tooth width
wes = awes*R;     % stator slot width
wa = wes;         % stator magnet width
Lper = (wa+wes+2*wds)*6; % magnetic period length equivalent to 180 degrees
ha  = 35e-3;      % magnet height
e   = 1.5e-3;     % Airgap
hs  = 35e-3;      % stator height
hes = 25e-3;      % slot height

hdr = 35e-03;     % height of rotor tooth
her = 25e-03;     % height of rotor slot 
wdr = awdr*R;     % width of rotor tooth
wer = awer*R;     % width of rotor slot

la  = 53e-03;     % active length
Br0 = 1.2;        % PM remanence
mu0 = 4*pi*1e-7;
mur1 = 1;
mur2 = 2000;      % Fer relative permeability


mds = 35;
mes = 25;
ma = mes;
mal = 13;
mar = 12;
p0s = 10;   
m  =  (mds*2+mes+ma)*6; 

p  = 14; 


% permeance unit dims
h  = hs/p;
l  = Lper/m;

% For the elements in the slots + airgap
Rv1 = (1/mur1/mu0)*l/(h*la);
Rw1 = (1/mur1/mu0)*h/(l*la);
Pv1 = mur1*mu0*(h*la/l);
Pw1 = mur1*mu0*(l*la/h);

% For the elements of the yoke and teeth
Rv20 = (1/mur2/mu0)*l/(h*la);
Rw20 = (1/mur2/mu0)*h/(l*la);
Pv20 = mur2*mu0*(h*la/l);
Pw20 = mur2*mu0*(l*la/h);

% for airgap fine meshing, 2 nodes in 1.5 mm
h1  = 1.5e-03/2; 
Rwag  = (1/mur1/mu0)*h1/(l*la);
Rvag = (1/mur1/mu0)*l/(h1*la);


% Rotor Horizontal elements
m0m = 30;
mdr = m0m*2; % m0 rotor tooth
mer = 84; % m0 rotor slot
mr = 5*(mdr+mer); % mr == m
p0m = 1;
pr  = (14)*p0m;
p1 = pr+2;

pcr = 4*p0m;

nn = (p+p1)*m;     % Total number of nodes of the RN
global M nnum ;
M = zeros(nn);
MP = zeros(nn);
nnum = zeros((p+p1),m);

%% rotor culasse meshing 
% First line
r11 = 1;
s11 = 1;
nnum(r11,s11)=(r11-1)*mr+s11;
MP(nnum(r11,s11),nnum(r11,s11))     =  1/(Rw20/2);
MP(nnum(r11,s11),nnum(r11,s11)+m-1) =  1/Rv20; 
MP(nnum(r11,s11),nnum(r11,s11)+1)   =  1/Rv20;
MP(nnum(r11,s11),nnum(r11,s11)+m)   =  1/Rw20;
M(nnum(r11,s11),nnum(r11,s11))      = MP(nnum(r11,s11),nnum(r11,s11))+MP(nnum(r11,s11),nnum(r11,s11)+m-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m);
M(nnum(r11,s11),nnum(r11,s11)+m-1)  = -MP(nnum(r11,s11),nnum(r11,s11)+m-1);
M(nnum(r11,s11),nnum(r11,s11)+1)    = -MP(nnum(r11,s11),nnum(r11,s11)+1);
M(nnum(r11,s11),nnum(r11,s11)+m)    = -MP(nnum(r11,s11),nnum(r11,s11)+m);
for s11 = 2:mr-1
    nnum(r11,s11)=(r11-1)*mr+s11;
    MP(nnum(r11,s11),nnum(r11,s11))   = 1/(Rw20/2);
    MP(nnum(r11,s11),nnum(r11,s11)-1) = 1/Rv20;
    MP(nnum(r11,s11),nnum(r11,s11)+1) = 1/Rv20;
    MP(nnum(r11,s11),nnum(r11,s11)+m) = 1/Rw20;
    M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11))+MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m);
    M(nnum(r11,s11),nnum(r11,s11)-1)  = -MP(nnum(r11,s11),nnum(r11,s11)-1);
    M(nnum(r11,s11),nnum(r11,s11)+1)  = -MP(nnum(r11,s11),nnum(r11,s11)+1);
    M(nnum(r11,s11),nnum(r11,s11)+m)  = -MP(nnum(r11,s11),nnum(r11,s11)+m);
end
s11 = mr;
nnum(r11,s11)=(r11-1)*mr+s11;
MP(nnum(r11,s11),nnum(r11,s11))     =  1/(Rw20/2);
MP(nnum(r11,s11),nnum(r11,s11)-1)   =  1/Rv20;
MP(nnum(r11,s11),nnum(r11,s11)+m)   =  1/Rw20;
MP(nnum(r11,s11),nnum(r11,s11)-m+1) =  1/Rv20;
M(nnum(r11,s11),nnum(r11,s11))      = MP(nnum(r11,s11),nnum(r11,s11))+MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m+1);
M(nnum(r11,s11),nnum(r11,s11)-1)    =  -MP(nnum(r11,s11),nnum(r11,s11)-1);
M(nnum(r11,s11),nnum(r11,s11)+m)    =  -MP(nnum(r11,s11),nnum(r11,s11)+m);
M(nnum(r11,s11),nnum(r11,s11)-m+1)  =  -MP(nnum(r11,s11),nnum(r11,s11)-m+1);

% first coulumn of rotor culasse
s11 = 1;
for r11 = 2:pcr-1
    nnum(r11,s11)=(r11-1)*mr+s11;
    MP(nnum(r11,s11),nnum(r11,s11)+m-1)   =   1/Rv20;
    MP(nnum(r11,s11),nnum(r11,s11)+1)     =   1/Rv20;
    MP(nnum(r11,s11),nnum(r11,s11)+m)     =   1/Rw20;
    MP(nnum(r11,s11),nnum(r11,s11)-m)     =   1/Rw20;
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
    MP(nnum(r11,s11),nnum(r11,s11)-m+1)   =  1/Rv20;
    MP(nnum(r11,s11),nnum(r11,s11)-1)     =   1/Rv20;
    MP(nnum(r11,s11),nnum(r11,s11)+m)     =   1/Rw20;
    MP(nnum(r11,s11),nnum(r11,s11)-m)     =   1/Rw20;
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
        MP(nnum(r11,s11),nnum(r11,s11)-1) = 1/Rv20;
        MP(nnum(r11,s11),nnum(r11,s11)+1) = 1/Rv20;
        MP(nnum(r11,s11),nnum(r11,s11)+m) = 1/Rw20;
        MP(nnum(r11,s11),nnum(r11,s11)-m) = 1/Rw20;
        M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
        M(nnum(r11,s11),nnum(r11,s11)-1) = -MP(nnum(r11,s11),nnum(r11,s11)-1);
        M(nnum(r11,s11),nnum(r11,s11)+1) = -MP(nnum(r11,s11),nnum(r11,s11)+1);
        M(nnum(r11,s11),nnum(r11,s11)+m) = -MP(nnum(r11,s11),nnum(r11,s11)+m);
        M(nnum(r11,s11),nnum(r11,s11)-m) = -MP(nnum(r11,s11),nnum(r11,s11)-m);
    end
end

%% Border culasse rotor/teeth rotor
% row number pcr
r11 = pcr;
s11 = 1;
nnum(r11,s11)=(r11-1)*mr+s11;
MP(nnum(r11,s11),nnum(r11,s11)-m)   =  1/Rw20;
MP(nnum(r11,s11),nnum(r11,s11)+m-1) =  1/Rv20;
MP(nnum(r11,s11),nnum(r11,s11)+1)   =  1/Rv20;
MP(nnum(r11,s11),nnum(r11,s11)+m)   =  1/Rw20; % permeance with the rotor tooth above
M(nnum(r11,s11),nnum(r11,s11))      = MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m-1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
M(nnum(r11,s11),nnum(r11,s11)-m)    =  -MP(nnum(r11,s11),nnum(r11,s11)-m);
M(nnum(r11,s11),nnum(r11,s11)+m-1)  =  -MP(nnum(r11,s11),nnum(r11,s11)+m-1);
M(nnum(r11,s11),nnum(r11,s11)+1)    =  -MP(nnum(r11,s11),nnum(r11,s11)+1);
M(nnum(r11,s11),nnum(r11,s11)+m)    =  -MP(nnum(r11,s11),nnum(r11,s11)+m);
for s11 = 2:m0m
    nnum(r11,s11)=(r11-1)*mr+s11;
    MP(nnum(r11,s11),nnum(r11,s11)-m) = 1/Rw20;
    MP(nnum(r11,s11),nnum(r11,s11)-1) = 1/Rv20;
    MP(nnum(r11,s11),nnum(r11,s11)+1) = 1/Rv20;
    MP(nnum(r11,s11),nnum(r11,s11)+m) = 1/(Rw20);
    M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
    M(nnum(r11,s11),nnum(r11,s11)-1) = -MP(nnum(r11,s11),nnum(r11,s11)-1);
    M(nnum(r11,s11),nnum(r11,s11)+1) = -MP(nnum(r11,s11),nnum(r11,s11)+1);
    M(nnum(r11,s11),nnum(r11,s11)+m) = -MP(nnum(r11,s11),nnum(r11,s11)+m);
    M(nnum(r11,s11),nnum(r11,s11)-m) = -MP(nnum(r11,s11),nnum(r11,s11)-m);
end
for dent = 1:4
    for s11 = (mer+m0m)+1+(dent-1)*(mer+mdr):(mer+m0m+mdr)+(dent-1)*(mer+mdr)
        nnum(r11,s11)=(r11-1)*m+s11;
        MP(nnum(r11,s11),nnum(r11,s11)-m) = 1/Rw20;
        MP(nnum(r11,s11),nnum(r11,s11)-1) = 1/Rv20;
        MP(nnum(r11,s11),nnum(r11,s11)+1) = 1/Rv20;
        MP(nnum(r11,s11),nnum(r11,s11)+m) = 1/Rw20;
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
        MP(nnum(r11,s11),nnum(r11,s11)-m) = 1/Rw20;
        MP(nnum(r11,s11),nnum(r11,s11)-1) = 1/Rv20;
        MP(nnum(r11,s11),nnum(r11,s11)+1) = 1/Rv20;
        MP(nnum(r11,s11),nnum(r11,s11)+m) = 1/(Rw1/2+Rw20/2); % slot air above
        M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
        M(nnum(r11,s11),nnum(r11,s11)-1) = -MP(nnum(r11,s11),nnum(r11,s11)-1);
        M(nnum(r11,s11),nnum(r11,s11)+1) = -MP(nnum(r11,s11),nnum(r11,s11)+1);
        M(nnum(r11,s11),nnum(r11,s11)+m) = -MP(nnum(r11,s11),nnum(r11,s11)+m);
        M(nnum(r11,s11),nnum(r11,s11)-m) = -MP(nnum(r11,s11),nnum(r11,s11)-m);
    end
end
for s11 = (5*mer+4*mdr+m0m)+1:mr-1 % here we reached the last half stator tooth
    nnum(r11,s11)=(r11-1)*mr+s11;
    MP(nnum(r11,s11),nnum(r11,s11)-m) = 1/Rw20;
    MP(nnum(r11,s11),nnum(r11,s11)-1) = 1/Rv20;
    MP(nnum(r11,s11),nnum(r11,s11)+1) = 1/Rv20;
    MP(nnum(r11,s11),nnum(r11,s11)+m) = 1/Rw20;
    M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
    M(nnum(r11,s11),nnum(r11,s11)-1) = -MP(nnum(r11,s11),nnum(r11,s11)-1);
    M(nnum(r11,s11),nnum(r11,s11)+1) = -MP(nnum(r11,s11),nnum(r11,s11)+1);
    M(nnum(r11,s11),nnum(r11,s11)+m) = -MP(nnum(r11,s11),nnum(r11,s11)+m);
    M(nnum(r11,s11),nnum(r11,s11)-m) = -MP(nnum(r11,s11),nnum(r11,s11)-m);
end
s11 = mr;
nnum(r11,s11)=(r11-1)*mr+s11;
MP(nnum(r11,s11),nnum(r11,s11)-m)   =  1/Rw20;
MP(nnum(r11,s11),nnum(r11,s11)-1)   =  1/Rv20;
MP(nnum(r11,s11),nnum(r11,s11)+m)   =  1/(Rw20);
MP(nnum(r11,s11),nnum(r11,s11)-m+1) =  1/Rv20;
M(nnum(r11,s11),nnum(r11,s11))      = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)-m+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
M(nnum(r11,s11),nnum(r11,s11)-1)   = -MP(nnum(r11,s11),nnum(r11,s11)-1);
M(nnum(r11,s11),nnum(r11,s11)-m+1) = -MP(nnum(r11,s11),nnum(r11,s11)-m+1);
M(nnum(r11,s11),nnum(r11,s11)+m)   = -MP(nnum(r11,s11),nnum(r11,s11)+m);
M(nnum(r11,s11),nnum(r11,s11)-m)   = -MP(nnum(r11,s11),nnum(r11,s11)-m);

%% rotor teeth and slots meshing r11 >= pcr+1
% first column s11 = 1
s11 = 1;
for r11 = pcr+1:pr-1
    nnum(r11,s11)=(r11-1)*m+s11;
    MP(nnum(r11,s11),nnum(r11,s11)+m-1) =  1/Rv20;
    MP(nnum(r11,s11),nnum(r11,s11)+1)   =  1/Rv20;
    MP(nnum(r11,s11),nnum(r11,s11)-m)   =  1/Rw20;
    MP(nnum(r11,s11),nnum(r11,s11)+m)   =  1/Rw20;
    M(nnum(r11,s11),nnum(r11,s11))      =  MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m-1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
    M(nnum(r11,s11),nnum(r11,s11)+m-1)  = -MP(nnum(r11,s11),nnum(r11,s11)+m-1);
    M(nnum(r11,s11),nnum(r11,s11)+1)    = -MP(nnum(r11,s11),nnum(r11,s11)+1);
    M(nnum(r11,s11),nnum(r11,s11)+m)    = -MP(nnum(r11,s11),nnum(r11,s11)+m);
    M(nnum(r11,s11),nnum(r11,s11)-m)    = -MP(nnum(r11,s11),nnum(r11,s11)-m);
end
% last row of rotor, airgap border
r11 = pr;
nnum(r11,s11)=(r11-1)*mr+s11;
MP(nnum(r11,s11),nnum(r11,s11)+m-1) =  1/Rv20;
MP(nnum(r11,s11),nnum(r11,s11)+1)   =  1/Rv20;
MP(nnum(r11,s11),nnum(r11,s11)-m)   =  1/Rw20;
MP(nnum(r11,s11),nnum(r11,s11)+m)   =  1/(Rw20/2 + Rwag/2);
M(nnum(r11,s11),nnum(r11,s11))      =  MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m-1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
M(nnum(r11,s11),nnum(r11,s11)+m-1)  = -MP(nnum(r11,s11),nnum(r11,s11)+m-1);
M(nnum(r11,s11),nnum(r11,s11)+1)    = -MP(nnum(r11,s11),nnum(r11,s11)+1);
M(nnum(r11,s11),nnum(r11,s11)+m)    = -MP(nnum(r11,s11),nnum(r11,s11)+m);
M(nnum(r11,s11),nnum(r11,s11)-m)    = -MP(nnum(r11,s11),nnum(r11,s11)-m);

r11 = pr+1; % first row of airgap
nnum(r11,s11)=(r11-1)*mr+s11;
MP(nnum(r11,s11),nnum(r11,s11)+m-1) =  1/Rvag;
MP(nnum(r11,s11),nnum(r11,s11)+1)   =  1/Rvag;
MP(nnum(r11,s11),nnum(r11,s11)+m)   =  1/Rwag;
MP(nnum(r11,s11),nnum(r11,s11)-m)   =  1/(Rw20/2 + Rwag/2);
M(nnum(r11,s11),nnum(r11,s11))      =  MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m-1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
M(nnum(r11,s11),nnum(r11,s11)+m-1)  = -MP(nnum(r11,s11),nnum(r11,s11)+m-1);
M(nnum(r11,s11),nnum(r11,s11)+1)    = -MP(nnum(r11,s11),nnum(r11,s11)+1);
M(nnum(r11,s11),nnum(r11,s11)+m)    = -MP(nnum(r11,s11),nnum(r11,s11)+m);
M(nnum(r11,s11),nnum(r11,s11)-m)    = -MP(nnum(r11,s11),nnum(r11,s11)-m);

% Last column of rotor (s11 = mr), r11 >= pcr+1
s11 = mr;
for r11 = pcr+1:pr-1
    nnum(r11,s11)=(r11-1)*mr+s11;
    MP(nnum(r11,s11),nnum(r11,s11)-m+1) =  1/Rv20;
    MP(nnum(r11,s11),nnum(r11,s11)-1)   =  1/Rv20;
    MP(nnum(r11,s11),nnum(r11,s11)-m)   =  1/Rw20;
    MP(nnum(r11,s11),nnum(r11,s11)+m)   =  1/Rw20;
    M(nnum(r11,s11),nnum(r11,s11))      = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)-m+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
    M(nnum(r11,s11),nnum(r11,s11)-m+1)  = -MP(nnum(r11,s11),nnum(r11,s11)-m+1);
    M(nnum(r11,s11),nnum(r11,s11)-1)    = -MP(nnum(r11,s11),nnum(r11,s11)-1);
    M(nnum(r11,s11),nnum(r11,s11)+m)    = -MP(nnum(r11,s11),nnum(r11,s11)+m);
    M(nnum(r11,s11),nnum(r11,s11)-m)    = -MP(nnum(r11,s11),nnum(r11,s11)-m);
end
r11 = pr;
nnum(r11,s11)=(r11-1)*mr+s11;
MP(nnum(r11,s11),nnum(r11,s11)-m+1) =  1/Rv20;
MP(nnum(r11,s11),nnum(r11,s11)-1)   =  1/Rv20;
MP(nnum(r11,s11),nnum(r11,s11)-m)   =  1/Rw20;
MP(nnum(r11,s11),nnum(r11,s11)+m)   =  1/(Rw20/2 + Rwag/2);
M(nnum(r11,s11),nnum(r11,s11))      = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)-m+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
M(nnum(r11,s11),nnum(r11,s11)-m+1)  = -MP(nnum(r11,s11),nnum(r11,s11)-m+1);
M(nnum(r11,s11),nnum(r11,s11)-1)    = -MP(nnum(r11,s11),nnum(r11,s11)-1);
M(nnum(r11,s11),nnum(r11,s11)+m)    = -MP(nnum(r11,s11),nnum(r11,s11)+m);
M(nnum(r11,s11),nnum(r11,s11)-m)    = -MP(nnum(r11,s11),nnum(r11,s11)-m);

r11 = pr+1;
nnum(r11,s11)=(r11-1)*mr+s11;
MP(nnum(r11,s11),nnum(r11,s11)-m+1) =  1/Rvag;
MP(nnum(r11,s11),nnum(r11,s11)-1)   =  1/Rvag;
MP(nnum(r11,s11),nnum(r11,s11)+m)   =  1/Rwag;
MP(nnum(r11,s11),nnum(r11,s11)-m)   =  1/(Rw20/2 + Rwag/2);
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
        MP(nnum(r11,s11),nnum(r11,s11)-1) = 1/Rv20;
        MP(nnum(r11,s11),nnum(r11,s11)+1) = 1/(Rv20/2+Rv1/2);
        MP(nnum(r11,s11),nnum(r11,s11)-m) = 1/Rw20;
        MP(nnum(r11,s11),nnum(r11,s11)+m) = 1/Rw20;
        M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
        M(nnum(r11,s11),nnum(r11,s11)-1)    = -MP(nnum(r11,s11),nnum(r11,s11)-1);
        M(nnum(r11,s11),nnum(r11,s11)+1)    = -MP(nnum(r11,s11),nnum(r11,s11)+1);
        M(nnum(r11,s11),nnum(r11,s11)+m)    = -MP(nnum(r11,s11),nnum(r11,s11)+m);
        M(nnum(r11,s11),nnum(r11,s11)-m)    = -MP(nnum(r11,s11),nnum(r11,s11)-m);
    end
    r11 = pr;
    nnum(r11,s11)=(r11-1)*mr+s11;
    MP(nnum(r11,s11),nnum(r11,s11)+1)   =  1/(Rv20/2+Rv1/2);
    MP(nnum(r11,s11),nnum(r11,s11)-1)   =  1/Rv20;
    MP(nnum(r11,s11),nnum(r11,s11)-m)   =  1/Rw20;
    MP(nnum(r11,s11),nnum(r11,s11)+m)   =  1/(Rw20/2 + Rwag/2);
    M(nnum(r11,s11),nnum(r11,s11))      = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
    M(nnum(r11,s11),nnum(r11,s11)+1)    = -MP(nnum(r11,s11),nnum(r11,s11)+1);
    M(nnum(r11,s11),nnum(r11,s11)-1)    = -MP(nnum(r11,s11),nnum(r11,s11)-1);
    M(nnum(r11,s11),nnum(r11,s11)+m)    = -MP(nnum(r11,s11),nnum(r11,s11)+m);
    M(nnum(r11,s11),nnum(r11,s11)-m)    = -MP(nnum(r11,s11),nnum(r11,s11)-m);
    r11 = pr+1;
    nnum(r11,s11)=(r11-1)*mr+s11;
    MP(nnum(r11,s11),nnum(r11,s11)+1)   =  1/Rvag;
    MP(nnum(r11,s11),nnum(r11,s11)-1)   =  1/Rvag;
    MP(nnum(r11,s11),nnum(r11,s11)+m)   =  1/Rwag;
    MP(nnum(r11,s11),nnum(r11,s11)-m)   =  1/(Rw20/2 + Rwag/2);
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
        MP(nnum(r11,s11),nnum(r11,s11)-1) = 1/(Rv20/2+Rv1/2);
        MP(nnum(r11,s11),nnum(r11,s11)+1) = 1/Rv20;
        MP(nnum(r11,s11),nnum(r11,s11)-m) = 1/Rw20;
        MP(nnum(r11,s11),nnum(r11,s11)+m) = 1/Rw20;
        M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
        M(nnum(r11,s11),nnum(r11,s11)-1)  = -MP(nnum(r11,s11),nnum(r11,s11)-1);
        M(nnum(r11,s11),nnum(r11,s11)+1)  = -MP(nnum(r11,s11),nnum(r11,s11)+1);
        M(nnum(r11,s11),nnum(r11,s11)+m)  = -MP(nnum(r11,s11),nnum(r11,s11)+m);
        M(nnum(r11,s11),nnum(r11,s11)-m)  = -MP(nnum(r11,s11),nnum(r11,s11)-m);
    end
    r11 = pr;
    nnum(r11,s11)=(r11-1)*mr+s11;
    MP(nnum(r11,s11),nnum(r11,s11)-1)   =  1/(Rv20/2+Rv1/2);
    MP(nnum(r11,s11),nnum(r11,s11)+1)   =  1/Rv20;
    MP(nnum(r11,s11),nnum(r11,s11)-m)   =  1/Rw20;
    MP(nnum(r11,s11),nnum(r11,s11)+m)   =  1/(Rw20/2 + Rwag/2);
    M(nnum(r11,s11),nnum(r11,s11))      = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
    M(nnum(r11,s11),nnum(r11,s11)+1)    = -MP(nnum(r11,s11),nnum(r11,s11)+1);
    M(nnum(r11,s11),nnum(r11,s11)-1)    = -MP(nnum(r11,s11),nnum(r11,s11)-1);
    M(nnum(r11,s11),nnum(r11,s11)+m)    = -MP(nnum(r11,s11),nnum(r11,s11)+m);
    M(nnum(r11,s11),nnum(r11,s11)-m)    = -MP(nnum(r11,s11),nnum(r11,s11)-m);
    r11 = pr+1;
    nnum(r11,s11)=(r11-1)*mr+s11;
    MP(nnum(r11,s11),nnum(r11,s11)-1)   =  1/Rvag;
    MP(nnum(r11,s11),nnum(r11,s11)+1)   =  1/Rvag;
    MP(nnum(r11,s11),nnum(r11,s11)+m)   =  1/Rwag;
    MP(nnum(r11,s11),nnum(r11,s11)-m)   =  1/(Rw20/2 + Rwag/2);
    M(nnum(r11,s11),nnum(r11,s11))      = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
    M(nnum(r11,s11),nnum(r11,s11)+1)    = -MP(nnum(r11,s11),nnum(r11,s11)+1);
    M(nnum(r11,s11),nnum(r11,s11)-1)    = -MP(nnum(r11,s11),nnum(r11,s11)-1);
    M(nnum(r11,s11),nnum(r11,s11)+m)    = -MP(nnum(r11,s11),nnum(r11,s11)+m);
    M(nnum(r11,s11),nnum(r11,s11)-m)    = -MP(nnum(r11,s11),nnum(r11,s11)-m);
end

%% the Half teeth insides without boders
% needed in case m0m > 2
for dent = 1:2
    for s11 = 2+(dent-1)*(5*mer+4*mdr+m0m):m0m-1+(dent-1)*(5*mer+4*mdr+m0m)
        for r11=pcr+1:pr-1
            nnum(r11,s11)=(r11-1)*m+s11;
            MP(nnum(r11,s11),nnum(r11,s11)-1) = 1/Rv20;
            MP(nnum(r11,s11),nnum(r11,s11)+1) = 1/Rv20;
            MP(nnum(r11,s11),nnum(r11,s11)-m) = 1/Rw20;
            MP(nnum(r11,s11),nnum(r11,s11)+m) = 1/Rw20;
            M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
            M(nnum(r11,s11),nnum(r11,s11)-1)  = -MP(nnum(r11,s11),nnum(r11,s11)-1);
            M(nnum(r11,s11),nnum(r11,s11)+1)  = -MP(nnum(r11,s11),nnum(r11,s11)+1);
            M(nnum(r11,s11),nnum(r11,s11)+m)  = -MP(nnum(r11,s11),nnum(r11,s11)+m);
            M(nnum(r11,s11),nnum(r11,s11)-m)  = -MP(nnum(r11,s11),nnum(r11,s11)-m);
        end
    r11 = pr;
    nnum(r11,s11)=(r11-1)*m+s11;
    MP(nnum(r11,s11),nnum(r11,s11)-1)   =  1/(Rv20);
    MP(nnum(r11,s11),nnum(r11,s11)+1)   =  1/Rv20;
    MP(nnum(r11,s11),nnum(r11,s11)-m)   =  1/Rw20;
    MP(nnum(r11,s11),nnum(r11,s11)+m)   =  1/(Rw20/2 + Rwag/2);
    M(nnum(r11,s11),nnum(r11,s11))      = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
    M(nnum(r11,s11),nnum(r11,s11)+1)    = -MP(nnum(r11,s11),nnum(r11,s11)+1);
    M(nnum(r11,s11),nnum(r11,s11)-1)    = -MP(nnum(r11,s11),nnum(r11,s11)-1);
    M(nnum(r11,s11),nnum(r11,s11)+m)    = -MP(nnum(r11,s11),nnum(r11,s11)+m);
    M(nnum(r11,s11),nnum(r11,s11)-m)    = -MP(nnum(r11,s11),nnum(r11,s11)-m);
    r11 = pr+1;
    nnum(r11,s11)=(r11-1)*m+s11;
    MP(nnum(r11,s11),nnum(r11,s11)-1)   =  1/(Rvag);
    MP(nnum(r11,s11),nnum(r11,s11)+1)   =  1/Rvag;
    MP(nnum(r11,s11),nnum(r11,s11)+m)   =  1/Rwag;
    MP(nnum(r11,s11),nnum(r11,s11)-m)   =  1/(Rw20/2 + Rwag/2);
    M(nnum(r11,s11),nnum(r11,s11))      = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
    M(nnum(r11,s11),nnum(r11,s11)+1)    = -MP(nnum(r11,s11),nnum(r11,s11)+1);
    M(nnum(r11,s11),nnum(r11,s11)-1)    = -MP(nnum(r11,s11),nnum(r11,s11)-1);
    M(nnum(r11,s11),nnum(r11,s11)+m)    = -MP(nnum(r11,s11),nnum(r11,s11)+m);
    M(nnum(r11,s11),nnum(r11,s11)-m)    = -MP(nnum(r11,s11),nnum(r11,s11)-m);
    end
end

%% the inside nodes of the rotor teeth
for dent = 1:4
    for s11 = m0m+mer+2+(dent-1)*(mer+mdr):m0m+mer+mdr-1+(dent-1)*(mer+mdr)
        for r11=pcr+1:pr-1
            nnum(r11,s11)=(r11-1)*m+s11;
            MP(nnum(r11,s11),nnum(r11,s11)-1) = 1/Rv20;
            MP(nnum(r11,s11),nnum(r11,s11)+1) = 1/Rv20;
            MP(nnum(r11,s11),nnum(r11,s11)-m) = 1/Rw20;
            MP(nnum(r11,s11),nnum(r11,s11)+m) = 1/Rw20;
            M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
            M(nnum(r11,s11),nnum(r11,s11)-1)  = -MP(nnum(r11,s11),nnum(r11,s11)-1);
            M(nnum(r11,s11),nnum(r11,s11)+1)  = -MP(nnum(r11,s11),nnum(r11,s11)+1);
            M(nnum(r11,s11),nnum(r11,s11)+m)  = -MP(nnum(r11,s11),nnum(r11,s11)+m);
            M(nnum(r11,s11),nnum(r11,s11)-m)  = -MP(nnum(r11,s11),nnum(r11,s11)-m);
        end
     r11 = pr;
     nnum(r11,s11)=(r11-1)*m+s11;
     MP(nnum(r11,s11),nnum(r11,s11)-1)   =  1/(Rv20);
     MP(nnum(r11,s11),nnum(r11,s11)+1)   =  1/Rv20;
     MP(nnum(r11,s11),nnum(r11,s11)-m)   =  1/Rw20;
     MP(nnum(r11,s11),nnum(r11,s11)+m)   =  1/(Rw20/2 + Rwag/2);
     M(nnum(r11,s11),nnum(r11,s11))      = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
     M(nnum(r11,s11),nnum(r11,s11)+1)    = -MP(nnum(r11,s11),nnum(r11,s11)+1);
     M(nnum(r11,s11),nnum(r11,s11)-1)    = -MP(nnum(r11,s11),nnum(r11,s11)-1);
     M(nnum(r11,s11),nnum(r11,s11)+m)    = -MP(nnum(r11,s11),nnum(r11,s11)+m);
     M(nnum(r11,s11),nnum(r11,s11)-m)    = -MP(nnum(r11,s11),nnum(r11,s11)-m);
     r11 = pr+1;
     nnum(r11,s11)=(r11-1)*m+s11;
     MP(nnum(r11,s11),nnum(r11,s11)-1)   =  1/(Rvag);
     MP(nnum(r11,s11),nnum(r11,s11)+1)   =  1/Rvag;
     MP(nnum(r11,s11),nnum(r11,s11)+m)   =  1/Rwag;
     MP(nnum(r11,s11),nnum(r11,s11)-m)   =  1/(Rw20/2 + Rwag/2);
     M(nnum(r11,s11),nnum(r11,s11))      = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
     M(nnum(r11,s11),nnum(r11,s11)+1)    = -MP(nnum(r11,s11),nnum(r11,s11)+1);
     M(nnum(r11,s11),nnum(r11,s11)-1)    = -MP(nnum(r11,s11),nnum(r11,s11)-1);
     M(nnum(r11,s11),nnum(r11,s11)+m)    = -MP(nnum(r11,s11),nnum(r11,s11)+m);
     M(nnum(r11,s11),nnum(r11,s11)-m)    = -MP(nnum(r11,s11),nnum(r11,s11)-m);
    end
end

%% nodes inside slots in contact with the teeth
% left borders of the slots
for enc = 1:5
    s11 = m0m+1+(enc-1)*(mdr+mer);
    r11 = pcr+1;
    nnum(r11,s11)=(r11-1)*m+s11;
    MP(nnum(r11,s11),nnum(r11,s11)-1) = 1/(Rv20/2+Rv1/2);
    MP(nnum(r11,s11),nnum(r11,s11)+1) = 1/Rv1;
    MP(nnum(r11,s11),nnum(r11,s11)-m) = 1/(Rw1/2+Rw20/2);
    MP(nnum(r11,s11),nnum(r11,s11)+m) = 1/Rw1;
    M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
    M(nnum(r11,s11),nnum(r11,s11)-1)  = -MP(nnum(r11,s11),nnum(r11,s11)-1);
    M(nnum(r11,s11),nnum(r11,s11)+1)  = -MP(nnum(r11,s11),nnum(r11,s11)+1);
    M(nnum(r11,s11),nnum(r11,s11)+m)  = -MP(nnum(r11,s11),nnum(r11,s11)+m);
    M(nnum(r11,s11),nnum(r11,s11)-m)  = -MP(nnum(r11,s11),nnum(r11,s11)-m);
    for r11 = pcr+2:pr-1
        nnum(r11,s11)=(r11-1)*m+s11;
        MP(nnum(r11,s11),nnum(r11,s11)-1) = 1/(Rv20/2+Rv1/2);
        MP(nnum(r11,s11),nnum(r11,s11)+1) = 1/Rv1;
        MP(nnum(r11,s11),nnum(r11,s11)-m) = 1/Rw1;
        MP(nnum(r11,s11),nnum(r11,s11)+m) = 1/Rw1;
        M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
        M(nnum(r11,s11),nnum(r11,s11)-1)  = -MP(nnum(r11,s11),nnum(r11,s11)-1);
        M(nnum(r11,s11),nnum(r11,s11)+1)  = -MP(nnum(r11,s11),nnum(r11,s11)+1);
        M(nnum(r11,s11),nnum(r11,s11)+m)  = -MP(nnum(r11,s11),nnum(r11,s11)+m);
        M(nnum(r11,s11),nnum(r11,s11)-m)  = -MP(nnum(r11,s11),nnum(r11,s11)-m);
    end
    r11 = pr;
    nnum(r11,s11)=(r11-1)*m+s11;
    MP(nnum(r11,s11),nnum(r11,s11)-1) = 1/(Rv20/2+Rv1/2);
    MP(nnum(r11,s11),nnum(r11,s11)+1) = 1/Rv1;
    MP(nnum(r11,s11),nnum(r11,s11)-m) = 1/Rw1;
    MP(nnum(r11,s11),nnum(r11,s11)+m) = 1/(Rw1/2+Rwag/2);
    M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
    M(nnum(r11,s11),nnum(r11,s11)-1)  = -MP(nnum(r11,s11),nnum(r11,s11)-1);
    M(nnum(r11,s11),nnum(r11,s11)+1)  = -MP(nnum(r11,s11),nnum(r11,s11)+1);
    M(nnum(r11,s11),nnum(r11,s11)+m)  = -MP(nnum(r11,s11),nnum(r11,s11)+m);
    M(nnum(r11,s11),nnum(r11,s11)-m)  = -MP(nnum(r11,s11),nnum(r11,s11)-m);
    r11 = pr+1;
    nnum(r11,s11)=(r11-1)*m+s11;
    MP(nnum(r11,s11),nnum(r11,s11)-1) = 1/(Rvag);
    MP(nnum(r11,s11),nnum(r11,s11)+1) = 1/Rvag;
    MP(nnum(r11,s11),nnum(r11,s11)+m) = 1/Rwag;
    MP(nnum(r11,s11),nnum(r11,s11)-m) = 1/(Rw1/2+Rwag/2);
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
    MP(nnum(r11,s11),nnum(r11,s11)+1) = 1/(Rv20/2+Rv1/2);
    MP(nnum(r11,s11),nnum(r11,s11)-1) = 1/Rv1;
    MP(nnum(r11,s11),nnum(r11,s11)-m) = 1/(Rw1/2+Rw20/2);
    MP(nnum(r11,s11),nnum(r11,s11)+m) = 1/Rw1;
    M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
    M(nnum(r11,s11),nnum(r11,s11)-1)  = -MP(nnum(r11,s11),nnum(r11,s11)-1);
    M(nnum(r11,s11),nnum(r11,s11)+1)  = -MP(nnum(r11,s11),nnum(r11,s11)+1);
    M(nnum(r11,s11),nnum(r11,s11)+m)  = -MP(nnum(r11,s11),nnum(r11,s11)+m);
    M(nnum(r11,s11),nnum(r11,s11)-m)  = -MP(nnum(r11,s11),nnum(r11,s11)-m);
    for r11 = pcr+2:pr-1
        nnum(r11,s11)=(r11-1)*m+s11;
        MP(nnum(r11,s11),nnum(r11,s11)+1) = 1/(Rv20/2+Rv1/2);
        MP(nnum(r11,s11),nnum(r11,s11)-1) = 1/Rv1;
        MP(nnum(r11,s11),nnum(r11,s11)-m) = 1/Rw1;
        MP(nnum(r11,s11),nnum(r11,s11)+m) = 1/Rw1;
        M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
        M(nnum(r11,s11),nnum(r11,s11)-1)  = -MP(nnum(r11,s11),nnum(r11,s11)-1);
        M(nnum(r11,s11),nnum(r11,s11)+1)  = -MP(nnum(r11,s11),nnum(r11,s11)+1);
        M(nnum(r11,s11),nnum(r11,s11)+m)  = -MP(nnum(r11,s11),nnum(r11,s11)+m);
        M(nnum(r11,s11),nnum(r11,s11)-m)  = -MP(nnum(r11,s11),nnum(r11,s11)-m);
    end
    r11 = pr;
    nnum(r11,s11)=(r11-1)*m+s11;
    MP(nnum(r11,s11),nnum(r11,s11)+1) = 1/(Rv20/2+Rv1/2);
    MP(nnum(r11,s11),nnum(r11,s11)-1) = 1/Rv1;
    MP(nnum(r11,s11),nnum(r11,s11)-m) = 1/Rw1;
    MP(nnum(r11,s11),nnum(r11,s11)+m) = 1/(Rw1/2+Rwag/2);
    M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
    M(nnum(r11,s11),nnum(r11,s11)-1)  = -MP(nnum(r11,s11),nnum(r11,s11)-1);
    M(nnum(r11,s11),nnum(r11,s11)+1)  = -MP(nnum(r11,s11),nnum(r11,s11)+1);
    M(nnum(r11,s11),nnum(r11,s11)+m)  = -MP(nnum(r11,s11),nnum(r11,s11)+m);
    M(nnum(r11,s11),nnum(r11,s11)-m)  = -MP(nnum(r11,s11),nnum(r11,s11)-m);
    r11 = pr+1;
    nnum(r11,s11)=(r11-1)*m+s11;
    MP(nnum(r11,s11),nnum(r11,s11)+1) = 1/(Rvag);
    MP(nnum(r11,s11),nnum(r11,s11)-1) = 1/Rvag;
    MP(nnum(r11,s11),nnum(r11,s11)+m) = 1/Rwag;
    MP(nnum(r11,s11),nnum(r11,s11)-m) = 1/(Rw1/2+Rwag/2);
    M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
    M(nnum(r11,s11),nnum(r11,s11)-1)  = -MP(nnum(r11,s11),nnum(r11,s11)-1);
    M(nnum(r11,s11),nnum(r11,s11)+1)  = -MP(nnum(r11,s11),nnum(r11,s11)+1);
    M(nnum(r11,s11),nnum(r11,s11)+m)  = -MP(nnum(r11,s11),nnum(r11,s11)+m);
    M(nnum(r11,s11),nnum(r11,s11)-m)  = -MP(nnum(r11,s11),nnum(r11,s11)-m);
    
end

%% inner nodes of the rotor slots
for enc = 1:5
    for s11 = m0m+2+(enc-1)*(mer+mdr):m0m+mer-1+(enc-1)*(mer+mdr)
        r11 = pcr+1;
        nnum(r11,s11)=(r11-1)*m+s11;
        MP(nnum(r11,s11),nnum(r11,s11)-1) = 1/Rv1;
        MP(nnum(r11,s11),nnum(r11,s11)+1) = 1/Rv1;
        MP(nnum(r11,s11),nnum(r11,s11)-m) = 1/(Rw1/2+Rw20/2);
        MP(nnum(r11,s11),nnum(r11,s11)+m) = 1/Rw1;
        M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
        M(nnum(r11,s11),nnum(r11,s11)-1)  = -MP(nnum(r11,s11),nnum(r11,s11)-1);
        M(nnum(r11,s11),nnum(r11,s11)+1)  = -MP(nnum(r11,s11),nnum(r11,s11)+1);
        M(nnum(r11,s11),nnum(r11,s11)+m)  = -MP(nnum(r11,s11),nnum(r11,s11)+m);
        M(nnum(r11,s11),nnum(r11,s11)-m)  = -MP(nnum(r11,s11),nnum(r11,s11)-m);
        for r11 = pcr+2:pr-1
            nnum(r11,s11)=(r11-1)*m+s11;
            MP(nnum(r11,s11),nnum(r11,s11)-1) = 1/Rv1;
            MP(nnum(r11,s11),nnum(r11,s11)+1) = 1/Rv1;
            MP(nnum(r11,s11),nnum(r11,s11)-m) = 1/Rw1;
            MP(nnum(r11,s11),nnum(r11,s11)+m) = 1/Rw1;
            M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
            M(nnum(r11,s11),nnum(r11,s11)-1)  = -MP(nnum(r11,s11),nnum(r11,s11)-1);
            M(nnum(r11,s11),nnum(r11,s11)+1)  = -MP(nnum(r11,s11),nnum(r11,s11)+1);
            M(nnum(r11,s11),nnum(r11,s11)+m)  = -MP(nnum(r11,s11),nnum(r11,s11)+m);
            M(nnum(r11,s11),nnum(r11,s11)-m)  = -MP(nnum(r11,s11),nnum(r11,s11)-m);
        end
        r11 = pr;
        nnum(r11,s11)=(r11-1)*m+s11;
        MP(nnum(r11,s11),nnum(r11,s11)-1) = 1/Rv1;
        MP(nnum(r11,s11),nnum(r11,s11)+1) = 1/Rv1;
        MP(nnum(r11,s11),nnum(r11,s11)-m) = 1/Rw1;
        MP(nnum(r11,s11),nnum(r11,s11)+m) = 1/(Rw1/2+Rwag/2);
        M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
        M(nnum(r11,s11),nnum(r11,s11)-1)  = -MP(nnum(r11,s11),nnum(r11,s11)-1);
        M(nnum(r11,s11),nnum(r11,s11)+1)  = -MP(nnum(r11,s11),nnum(r11,s11)+1);
        M(nnum(r11,s11),nnum(r11,s11)+m)  = -MP(nnum(r11,s11),nnum(r11,s11)+m);
        M(nnum(r11,s11),nnum(r11,s11)-m)  = -MP(nnum(r11,s11),nnum(r11,s11)-m);
        r11 = pr+1;
        nnum(r11,s11)=(r11-1)*m+s11;
        MP(nnum(r11,s11),nnum(r11,s11)-1) = 1/Rvag;
        MP(nnum(r11,s11),nnum(r11,s11)+1) = 1/Rvag;
        MP(nnum(r11,s11),nnum(r11,s11)+m) = 1/Rwag;
        MP(nnum(r11,s11),nnum(r11,s11)-m) = 1/(Rw1/2+Rwag/2);
        M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
        M(nnum(r11,s11),nnum(r11,s11)-1)  = -MP(nnum(r11,s11),nnum(r11,s11)-1);
        M(nnum(r11,s11),nnum(r11,s11)+1)  = -MP(nnum(r11,s11),nnum(r11,s11)+1);
        M(nnum(r11,s11),nnum(r11,s11)+m)  = -MP(nnum(r11,s11),nnum(r11,s11)+m);
        M(nnum(r11,s11),nnum(r11,s11)-m)  = -MP(nnum(r11,s11),nnum(r11,s11)-m);
    end
end

%% stator meshing
%% PM meshing

% first column s11 = 1, r11>=p1
s11 = 1;
r11 = p1;
nnum(r11,s11)=(r11-1)*m+s11;
MP(nnum(r11,s11),nnum(r11,s11)+m-1) =  1/Rvag;
MP(nnum(r11,s11),nnum(r11,s11)+1)   =  1/Rvag;
MP(nnum(r11,s11),nnum(r11,s11)+m)   =  1/(Rw1/2+Rwag/2);
MP(nnum(r11,s11),nnum(r11,s11)-m)   =  1/Rwag;
M(nnum(r11,s11),nnum(r11,s11))      = MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m-1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
M(nnum(r11,s11),nnum(r11,s11)+m-1)  = -MP(nnum(r11,s11),nnum(r11,s11)+m-1);
M(nnum(r11,s11),nnum(r11,s11)+1)    = -MP(nnum(r11,s11),nnum(r11,s11)+1);
M(nnum(r11,s11),nnum(r11,s11)+m)    = -MP(nnum(r11,s11),nnum(r11,s11)+m);
M(nnum(r11,s11),nnum(r11,s11)-m)    = -MP(nnum(r11,s11),nnum(r11,s11)-m);
r11 = p1+1;
nnum(r11,s11)=(r11-1)*m+s11;
MP(nnum(r11,s11),nnum(r11,s11)+m-1) =  1/Rv1;
MP(nnum(r11,s11),nnum(r11,s11)+1)   =  1/Rv1;
MP(nnum(r11,s11),nnum(r11,s11)-m)   =  1/(Rw1/2+Rwag/2);
MP(nnum(r11,s11),nnum(r11,s11)+m)   =  1/Rw1;
M(nnum(r11,s11),nnum(r11,s11))      = MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m-1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
M(nnum(r11,s11),nnum(r11,s11)+m-1)  = -MP(nnum(r11,s11),nnum(r11,s11)+m-1);
M(nnum(r11,s11),nnum(r11,s11)+1)    = -MP(nnum(r11,s11),nnum(r11,s11)+1);
M(nnum(r11,s11),nnum(r11,s11)+m)    = -MP(nnum(r11,s11),nnum(r11,s11)+m);
M(nnum(r11,s11),nnum(r11,s11)-m)    = -MP(nnum(r11,s11),nnum(r11,s11)-m);

for r11 = p1+2:p1+p-1
    nnum(r11,s11)=(r11-1)*m+s11;
    MP(nnum(r11,s11),nnum(r11,s11)+m-1) =  1/(Rv1);
    MP(nnum(r11,s11),nnum(r11,s11)+1)   =  1/Rv1;
    MP(nnum(r11,s11),nnum(r11,s11)-m)   =  1/Rw1;
    MP(nnum(r11,s11),nnum(r11,s11)+m)   =  1/Rw1;
    M(nnum(r11,s11),nnum(r11,s11))      = MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m-1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
    M(nnum(r11,s11),nnum(r11,s11)+m-1)  = -MP(nnum(r11,s11),nnum(r11,s11)+m-1);
    M(nnum(r11,s11),nnum(r11,s11)+1)    = -MP(nnum(r11,s11),nnum(r11,s11)+1);
    M(nnum(r11,s11),nnum(r11,s11)+m)    = -MP(nnum(r11,s11),nnum(r11,s11)+m);
    M(nnum(r11,s11),nnum(r11,s11)-m)    = -MP(nnum(r11,s11),nnum(r11,s11)-m);
end
r11 = p1+p;
nnum(r11,s11)=(r11-1)*m+s11;
MP(nnum(r11,s11),nnum(r11,s11)+m-1) =  1/(Rv1);
MP(nnum(r11,s11),nnum(r11,s11)+1)   =  1/Rv1;
MP(nnum(r11,s11),nnum(r11,s11)-m)   =  1/Rw1;
M(nnum(r11,s11),nnum(r11,s11))      = MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m-1)+MP(nnum(r11,s11),nnum(r11,s11)-m);
M(nnum(r11,s11),nnum(r11,s11)+m-1)  = -MP(nnum(r11,s11),nnum(r11,s11)+m-1);
M(nnum(r11,s11),nnum(r11,s11)+1)    = -MP(nnum(r11,s11),nnum(r11,s11)+1);
M(nnum(r11,s11),nnum(r11,s11)-m)    = -MP(nnum(r11,s11),nnum(r11,s11)-m);

% Last column (s11 = m), r11 >= p1 (magnet)
s11 = m;
r11 = p1;
nnum(r11,s11)=(r11-1)*m+s11;
MP(nnum(r11,s11),nnum(r11,s11)-m+1) =  1/Rvag;
MP(nnum(r11,s11),nnum(r11,s11)-1)   =  1/Rvag;
MP(nnum(r11,s11),nnum(r11,s11)+m)   =  1/(Rw1/2+Rwag/2);
MP(nnum(r11,s11),nnum(r11,s11)-m)   =  1/Rwag;
M(nnum(r11,s11),nnum(r11,s11))      = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)-m+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
M(nnum(r11,s11),nnum(r11,s11)-m+1)  = -MP(nnum(r11,s11),nnum(r11,s11)-m+1);
M(nnum(r11,s11),nnum(r11,s11)-1)    = -MP(nnum(r11,s11),nnum(r11,s11)-1);
M(nnum(r11,s11),nnum(r11,s11)+m)    = -MP(nnum(r11,s11),nnum(r11,s11)+m);
M(nnum(r11,s11),nnum(r11,s11)-m)    = -MP(nnum(r11,s11),nnum(r11,s11)-m);
r11 = p1+1;
nnum(r11,s11)=(r11-1)*m+s11;
MP(nnum(r11,s11),nnum(r11,s11)-m+1) =  1/(Rv1);
MP(nnum(r11,s11),nnum(r11,s11)-1)   =  1/Rv1;
MP(nnum(r11,s11),nnum(r11,s11)-m)   =  1/(Rw1/2+Rwag/2);
MP(nnum(r11,s11),nnum(r11,s11)+m)   =  1/Rw1;
M(nnum(r11,s11),nnum(r11,s11))      = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)-m+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
M(nnum(r11,s11),nnum(r11,s11)-m+1)  = -MP(nnum(r11,s11),nnum(r11,s11)-m+1);
M(nnum(r11,s11),nnum(r11,s11)-1)    = -MP(nnum(r11,s11),nnum(r11,s11)-1);
M(nnum(r11,s11),nnum(r11,s11)+m)    = -MP(nnum(r11,s11),nnum(r11,s11)+m);
M(nnum(r11,s11),nnum(r11,s11)-m)    = -MP(nnum(r11,s11),nnum(r11,s11)-m);
for r11 = p1+2:p1+p-1
    nnum(r11,s11)=(r11-1)*m+s11;
    MP(nnum(r11,s11),nnum(r11,s11)-m+1) =  1/(Rv1);
    MP(nnum(r11,s11),nnum(r11,s11)-1)   =  1/Rv1;
    MP(nnum(r11,s11),nnum(r11,s11)-m)   =  1/Rw1;
    MP(nnum(r11,s11),nnum(r11,s11)+m)   =  1/Rw1;
    M(nnum(r11,s11),nnum(r11,s11))      = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)-m+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
    M(nnum(r11,s11),nnum(r11,s11)-m+1)  = -MP(nnum(r11,s11),nnum(r11,s11)-m+1);
    M(nnum(r11,s11),nnum(r11,s11)-1)    = -MP(nnum(r11,s11),nnum(r11,s11)-1);
    M(nnum(r11,s11),nnum(r11,s11)+m)    = -MP(nnum(r11,s11),nnum(r11,s11)+m);
    M(nnum(r11,s11),nnum(r11,s11)-m)    = -MP(nnum(r11,s11),nnum(r11,s11)-m);
end
r11 = p1+p;
nnum(r11,s11)=(r11-1)*m+s11;
MP(nnum(r11,s11),nnum(r11,s11)-m+1) = +1/(Rv1);
MP(nnum(r11,s11),nnum(r11,s11)-1)   =  1/Rv1;
MP(nnum(r11,s11),nnum(r11,s11)-m)   =  1/Rw1;
M(nnum(r11,s11),nnum(r11,s11))      = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)-m+1)+MP(nnum(r11,s11),nnum(r11,s11)-m);
M(nnum(r11,s11),nnum(r11,s11)-m+1)  = -MP(nnum(r11,s11),nnum(r11,s11)-m+1);
M(nnum(r11,s11),nnum(r11,s11)-1)    = -MP(nnum(r11,s11),nnum(r11,s11)-1);
M(nnum(r11,s11),nnum(r11,s11)-m)    = -MP(nnum(r11,s11),nnum(r11,s11)-m);

% Columns right sides of magnets r11 >= p1
for dent = 1:6
    s11 = mal+(dent-1)*(2*mds+mes+ma);
    r11= p1;
    nnum(r11,s11)=(r11-1)*m+s11;
    MP(nnum(r11,s11),nnum(r11,s11)-1) = 1/Rvag;
    MP(nnum(r11,s11),nnum(r11,s11)+1) = 1/(Rvag);
    MP(nnum(r11,s11),nnum(r11,s11)+m) = 1/(Rw1/2+Rwag/2);
    MP(nnum(r11,s11),nnum(r11,s11)-m) = 1/Rwag;
    M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
    M(nnum(r11,s11),nnum(r11,s11)-1)    = -MP(nnum(r11,s11),nnum(r11,s11)-1);
    M(nnum(r11,s11),nnum(r11,s11)+1)    = -MP(nnum(r11,s11),nnum(r11,s11)+1);
    M(nnum(r11,s11),nnum(r11,s11)+m)    = -MP(nnum(r11,s11),nnum(r11,s11)+m);
    M(nnum(r11,s11),nnum(r11,s11)-m)    = -MP(nnum(r11,s11),nnum(r11,s11)-m);
    r11= p1+1;
    nnum(r11,s11)=(r11-1)*m+s11;
    MP(nnum(r11,s11),nnum(r11,s11)-1) = 1/Rv1;
    MP(nnum(r11,s11),nnum(r11,s11)+1) = 1/(Rv1/2+Rv20/2);
    MP(nnum(r11,s11),nnum(r11,s11)-m) = 1/(Rw1/2+Rwag/2);
    MP(nnum(r11,s11),nnum(r11,s11)+m) = 1/Rw1;
    M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
    M(nnum(r11,s11),nnum(r11,s11)-1)    = -MP(nnum(r11,s11),nnum(r11,s11)-1);
    M(nnum(r11,s11),nnum(r11,s11)+1)    = -MP(nnum(r11,s11),nnum(r11,s11)+1);
    M(nnum(r11,s11),nnum(r11,s11)+m)    = -MP(nnum(r11,s11),nnum(r11,s11)+m);
    M(nnum(r11,s11),nnum(r11,s11)-m)    = -MP(nnum(r11,s11),nnum(r11,s11)-m);
    for r11 = p1+2:p1+p-1
        nnum(r11,s11)=(r11-1)*m+s11;
        MP(nnum(r11,s11),nnum(r11,s11)-1) = 1/Rv1;
        MP(nnum(r11,s11),nnum(r11,s11)+1) = 1/(Rv20/2+Rv1/2);
        MP(nnum(r11,s11),nnum(r11,s11)-m) = 1/Rw1;
        MP(nnum(r11,s11),nnum(r11,s11)+m) = 1/Rw1;
        M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
        M(nnum(r11,s11),nnum(r11,s11)-1)    = -MP(nnum(r11,s11),nnum(r11,s11)-1);
        M(nnum(r11,s11),nnum(r11,s11)+1)    = -MP(nnum(r11,s11),nnum(r11,s11)+1);
        M(nnum(r11,s11),nnum(r11,s11)+m)    = -MP(nnum(r11,s11),nnum(r11,s11)+m);
        M(nnum(r11,s11),nnum(r11,s11)-m)    = -MP(nnum(r11,s11),nnum(r11,s11)-m);
    end
    r11 = p1+p;
    nnum(r11,s11)=(r11-1)*m+s11;
    MP(nnum(r11,s11),nnum(r11,s11)+1)   = +1/(Rv20/2+Rv1/2);
    MP(nnum(r11,s11),nnum(r11,s11)-1)   =  1/Rv1;
    MP(nnum(r11,s11),nnum(r11,s11)-m)   =  1/Rw1;
    M(nnum(r11,s11),nnum(r11,s11))      = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)-m);
    M(nnum(r11,s11),nnum(r11,s11)+1)    = -MP(nnum(r11,s11),nnum(r11,s11)+1);
    M(nnum(r11,s11),nnum(r11,s11)-1)    = -MP(nnum(r11,s11),nnum(r11,s11)-1);
    M(nnum(r11,s11),nnum(r11,s11)-m)    = -MP(nnum(r11,s11),nnum(r11,s11)-m);
end

% Columns left sides of magnets r11 >= p1
for dent = 1:6
    s11 = mal+(2*mds+mes)+1+(dent-1)*(2*mds+mes+ma);
    r11= p1;
    nnum(r11,s11)=(r11-1)*m+s11;
    MP(nnum(r11,s11),nnum(r11,s11)-1) = 1/Rvag;
    MP(nnum(r11,s11),nnum(r11,s11)+1) = 1/(Rvag);
    MP(nnum(r11,s11),nnum(r11,s11)+m) = 1/(Rw1/2+Rwag/2);
    MP(nnum(r11,s11),nnum(r11,s11)-m) = 1/Rwag;
    M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
    M(nnum(r11,s11),nnum(r11,s11)-1)    = -MP(nnum(r11,s11),nnum(r11,s11)-1);
    M(nnum(r11,s11),nnum(r11,s11)+1)    = -MP(nnum(r11,s11),nnum(r11,s11)+1);
    M(nnum(r11,s11),nnum(r11,s11)+m)    = -MP(nnum(r11,s11),nnum(r11,s11)+m);
    M(nnum(r11,s11),nnum(r11,s11)-m)    = -MP(nnum(r11,s11),nnum(r11,s11)-m);
    r11= p1+1;
    nnum(r11,s11)=(r11-1)*m+s11;
    MP(nnum(r11,s11),nnum(r11,s11)+1) = 1/Rv1;
    MP(nnum(r11,s11),nnum(r11,s11)-1) = 1/(Rv1/2+Rv20/2);
    MP(nnum(r11,s11),nnum(r11,s11)-m) = 1/(Rw1/2+Rwag/2);
    MP(nnum(r11,s11),nnum(r11,s11)+m) = 1/Rw1;
    M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
    M(nnum(r11,s11),nnum(r11,s11)-1)    = -MP(nnum(r11,s11),nnum(r11,s11)-1);
    M(nnum(r11,s11),nnum(r11,s11)+1)    = -MP(nnum(r11,s11),nnum(r11,s11)+1);
    M(nnum(r11,s11),nnum(r11,s11)+m)    = -MP(nnum(r11,s11),nnum(r11,s11)+m);
    M(nnum(r11,s11),nnum(r11,s11)-m)    = -MP(nnum(r11,s11),nnum(r11,s11)-m);
    for r11 = p1+2:p1+p-1
        nnum(r11,s11)=(r11-1)*m+s11;
        MP(nnum(r11,s11),nnum(r11,s11)+1) = 1/Rv1;
        MP(nnum(r11,s11),nnum(r11,s11)-1) = 1/(Rv20/2+Rv1/2);
        MP(nnum(r11,s11),nnum(r11,s11)-m) = 1/Rw1;
        MP(nnum(r11,s11),nnum(r11,s11)+m) = 1/Rw1;
        M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
        M(nnum(r11,s11),nnum(r11,s11)-1)    = -MP(nnum(r11,s11),nnum(r11,s11)-1);
        M(nnum(r11,s11),nnum(r11,s11)+1)    = -MP(nnum(r11,s11),nnum(r11,s11)+1);
        M(nnum(r11,s11),nnum(r11,s11)+m)    = -MP(nnum(r11,s11),nnum(r11,s11)+m);
        M(nnum(r11,s11),nnum(r11,s11)-m)    = -MP(nnum(r11,s11),nnum(r11,s11)-m);
    end
    r11 = p1+p;
    nnum(r11,s11)=(r11-1)*m+s11;
    MP(nnum(r11,s11),nnum(r11,s11)-1)   = +1/(Rv20/2+Rv1/2);
    MP(nnum(r11,s11),nnum(r11,s11)+1)   =  1/Rv1;
    MP(nnum(r11,s11),nnum(r11,s11)-m)   =  1/Rw1;
    M(nnum(r11,s11),nnum(r11,s11))      = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)-m);
    M(nnum(r11,s11),nnum(r11,s11)+1)    = -MP(nnum(r11,s11),nnum(r11,s11)+1);
    M(nnum(r11,s11),nnum(r11,s11)-1)    = -MP(nnum(r11,s11),nnum(r11,s11)-1);
    M(nnum(r11,s11),nnum(r11,s11)-m)    = -MP(nnum(r11,s11),nnum(r11,s11)-m);
end

%% the insides of the Half magnets
% needed in case m0s>2
for s11 = 2:mal-1
    r11= p1;
    nnum(r11,s11)=(r11-1)*m+s11;
    MP(nnum(r11,s11),nnum(r11,s11)-1) = 1/Rvag;
    MP(nnum(r11,s11),nnum(r11,s11)+1) = 1/(Rvag);
    MP(nnum(r11,s11),nnum(r11,s11)+m) = 1/(Rw1/2+Rwag/2);
    MP(nnum(r11,s11),nnum(r11,s11)-m) = 1/Rwag;
    M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
    M(nnum(r11,s11),nnum(r11,s11)-1)    = -MP(nnum(r11,s11),nnum(r11,s11)-1);
    M(nnum(r11,s11),nnum(r11,s11)+1)    = -MP(nnum(r11,s11),nnum(r11,s11)+1);
    M(nnum(r11,s11),nnum(r11,s11)+m)    = -MP(nnum(r11,s11),nnum(r11,s11)+m);
    M(nnum(r11,s11),nnum(r11,s11)-m)    = -MP(nnum(r11,s11),nnum(r11,s11)-m);
    r11= p1+1;
    nnum(r11,s11)=(r11-1)*m+s11;
    MP(nnum(r11,s11),nnum(r11,s11)-1) = 1/Rv1;
    MP(nnum(r11,s11),nnum(r11,s11)+1) = 1/(Rv1);
    MP(nnum(r11,s11),nnum(r11,s11)-m) = 1/(Rw1/2+Rwag/2);
    MP(nnum(r11,s11),nnum(r11,s11)+m) = 1/Rw1;
    M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
    M(nnum(r11,s11),nnum(r11,s11)-1)    = -MP(nnum(r11,s11),nnum(r11,s11)-1);
    M(nnum(r11,s11),nnum(r11,s11)+1)    = -MP(nnum(r11,s11),nnum(r11,s11)+1);
    M(nnum(r11,s11),nnum(r11,s11)+m)    = -MP(nnum(r11,s11),nnum(r11,s11)+m);
    M(nnum(r11,s11),nnum(r11,s11)-m)    = -MP(nnum(r11,s11),nnum(r11,s11)-m);
    for r11 = p1+2:p1+p-1
        nnum(r11,s11)=(r11-1)*m+s11;
        MP(nnum(r11,s11),nnum(r11,s11)+1) = 1/Rv1;
        MP(nnum(r11,s11),nnum(r11,s11)-1) = 1/Rv1;
        MP(nnum(r11,s11),nnum(r11,s11)-m) = 1/Rw1;
        MP(nnum(r11,s11),nnum(r11,s11)+m) = 1/Rw1;
        M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
        M(nnum(r11,s11),nnum(r11,s11)-1)    = -MP(nnum(r11,s11),nnum(r11,s11)-1);
        M(nnum(r11,s11),nnum(r11,s11)+1)    = -MP(nnum(r11,s11),nnum(r11,s11)+1);
        M(nnum(r11,s11),nnum(r11,s11)+m)    = -MP(nnum(r11,s11),nnum(r11,s11)+m);
        M(nnum(r11,s11),nnum(r11,s11)-m)    = -MP(nnum(r11,s11),nnum(r11,s11)-m);
    end
    r11 = p1+p;
    nnum(r11,s11)=(r11-1)*m+s11;
    MP(nnum(r11,s11),nnum(r11,s11)-1)   = +1/Rv1;
    MP(nnum(r11,s11),nnum(r11,s11)+1)   =  1/Rv1;
    MP(nnum(r11,s11),nnum(r11,s11)-m)   =  1/Rw1;
    M(nnum(r11,s11),nnum(r11,s11))      = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)-m);
    M(nnum(r11,s11),nnum(r11,s11)+1)    = -MP(nnum(r11,s11),nnum(r11,s11)+1);
    M(nnum(r11,s11),nnum(r11,s11)-1)    = -MP(nnum(r11,s11),nnum(r11,s11)-1);
    M(nnum(r11,s11),nnum(r11,s11)-m)    = -MP(nnum(r11,s11),nnum(r11,s11)-m);
end
for s11 = m-mar+2:m-1
    r11= p1;
    nnum(r11,s11)=(r11-1)*m+s11;
    MP(nnum(r11,s11),nnum(r11,s11)-1) = 1/Rvag;
    MP(nnum(r11,s11),nnum(r11,s11)+1) = 1/(Rvag);
    MP(nnum(r11,s11),nnum(r11,s11)+m) = 1/(Rw1/2+Rwag/2);
    MP(nnum(r11,s11),nnum(r11,s11)-m) = 1/Rwag;
    M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
    M(nnum(r11,s11),nnum(r11,s11)-1)    = -MP(nnum(r11,s11),nnum(r11,s11)-1);
    M(nnum(r11,s11),nnum(r11,s11)+1)    = -MP(nnum(r11,s11),nnum(r11,s11)+1);
    M(nnum(r11,s11),nnum(r11,s11)+m)    = -MP(nnum(r11,s11),nnum(r11,s11)+m);
    M(nnum(r11,s11),nnum(r11,s11)-m)    = -MP(nnum(r11,s11),nnum(r11,s11)-m);
    r11= p1+1;
    nnum(r11,s11)=(r11-1)*m+s11;
    MP(nnum(r11,s11),nnum(r11,s11)-1) = 1/Rv1;
    MP(nnum(r11,s11),nnum(r11,s11)+1) = 1/(Rv1);
    MP(nnum(r11,s11),nnum(r11,s11)-m) = 1/(Rw1/2+Rwag/2);
    MP(nnum(r11,s11),nnum(r11,s11)+m) = 1/Rw1;
    M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
    M(nnum(r11,s11),nnum(r11,s11)-1)    = -MP(nnum(r11,s11),nnum(r11,s11)-1);
    M(nnum(r11,s11),nnum(r11,s11)+1)    = -MP(nnum(r11,s11),nnum(r11,s11)+1);
    M(nnum(r11,s11),nnum(r11,s11)+m)    = -MP(nnum(r11,s11),nnum(r11,s11)+m);
    M(nnum(r11,s11),nnum(r11,s11)-m)    = -MP(nnum(r11,s11),nnum(r11,s11)-m);
    for r11 = p1+2:p1+p-1
        nnum(r11,s11)=(r11-1)*m+s11;
        MP(nnum(r11,s11),nnum(r11,s11)+1) = 1/Rv1;
        MP(nnum(r11,s11),nnum(r11,s11)-1) = 1/Rv1;
        MP(nnum(r11,s11),nnum(r11,s11)-m) = 1/Rw1;
        MP(nnum(r11,s11),nnum(r11,s11)+m) = 1/Rw1;
        M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
        M(nnum(r11,s11),nnum(r11,s11)-1)    = -MP(nnum(r11,s11),nnum(r11,s11)-1);
        M(nnum(r11,s11),nnum(r11,s11)+1)    = -MP(nnum(r11,s11),nnum(r11,s11)+1);
        M(nnum(r11,s11),nnum(r11,s11)+m)    = -MP(nnum(r11,s11),nnum(r11,s11)+m);
        M(nnum(r11,s11),nnum(r11,s11)-m)    = -MP(nnum(r11,s11),nnum(r11,s11)-m);
    end
    r11 = p1+p;
    nnum(r11,s11)=(r11-1)*m+s11;
    MP(nnum(r11,s11),nnum(r11,s11)-1)   = +1/Rv1;
    MP(nnum(r11,s11),nnum(r11,s11)+1)   =  1/Rv1;
    MP(nnum(r11,s11),nnum(r11,s11)-m)   =  1/Rw1;
    M(nnum(r11,s11),nnum(r11,s11))      = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)-m);
    M(nnum(r11,s11),nnum(r11,s11)+1)    = -MP(nnum(r11,s11),nnum(r11,s11)+1);
    M(nnum(r11,s11),nnum(r11,s11)-1)    = -MP(nnum(r11,s11),nnum(r11,s11)-1);
    M(nnum(r11,s11),nnum(r11,s11)-m)    = -MP(nnum(r11,s11),nnum(r11,s11)-m);
end

%% the inside nodes of the magnets
for dent = 1:5
    for s11 = mal+(2*mds+mes)+2+(dent-1)*(2*mds+mes+ma):mal+(2*mds+mes+ma)-1+(dent-1)*(2*mds+mes+ma)
    r11= p1;
    nnum(r11,s11)=(r11-1)*m+s11;
    MP(nnum(r11,s11),nnum(r11,s11)-1) = 1/Rvag;
    MP(nnum(r11,s11),nnum(r11,s11)+1) = 1/(Rvag);
    MP(nnum(r11,s11),nnum(r11,s11)+m) = 1/(Rw1/2+Rwag/2);
    MP(nnum(r11,s11),nnum(r11,s11)-m) = 1/Rwag;
    M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
    M(nnum(r11,s11),nnum(r11,s11)-1)    = -MP(nnum(r11,s11),nnum(r11,s11)-1);
    M(nnum(r11,s11),nnum(r11,s11)+1)    = -MP(nnum(r11,s11),nnum(r11,s11)+1);
    M(nnum(r11,s11),nnum(r11,s11)+m)    = -MP(nnum(r11,s11),nnum(r11,s11)+m);
    M(nnum(r11,s11),nnum(r11,s11)-m)    = -MP(nnum(r11,s11),nnum(r11,s11)-m);
    r11= p1+1;
    nnum(r11,s11)=(r11-1)*m+s11;
    MP(nnum(r11,s11),nnum(r11,s11)-1) = 1/Rv1;
    MP(nnum(r11,s11),nnum(r11,s11)+1) = 1/(Rv1);
    MP(nnum(r11,s11),nnum(r11,s11)-m) = 1/(Rw1/2+Rwag/2);
    MP(nnum(r11,s11),nnum(r11,s11)+m) = 1/Rw1;
    M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
    M(nnum(r11,s11),nnum(r11,s11)-1)    = -MP(nnum(r11,s11),nnum(r11,s11)-1);
    M(nnum(r11,s11),nnum(r11,s11)+1)    = -MP(nnum(r11,s11),nnum(r11,s11)+1);
    M(nnum(r11,s11),nnum(r11,s11)+m)    = -MP(nnum(r11,s11),nnum(r11,s11)+m);
    M(nnum(r11,s11),nnum(r11,s11)-m)    = -MP(nnum(r11,s11),nnum(r11,s11)-m);
    for r11 = p1+2:p1+p-1
        nnum(r11,s11)=(r11-1)*m+s11;
        MP(nnum(r11,s11),nnum(r11,s11)+1) = 1/Rv1;
        MP(nnum(r11,s11),nnum(r11,s11)-1) = 1/Rv1;
        MP(nnum(r11,s11),nnum(r11,s11)-m) = 1/Rw1;
        MP(nnum(r11,s11),nnum(r11,s11)+m) = 1/Rw1;
        M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
        M(nnum(r11,s11),nnum(r11,s11)-1)    = -MP(nnum(r11,s11),nnum(r11,s11)-1);
        M(nnum(r11,s11),nnum(r11,s11)+1)    = -MP(nnum(r11,s11),nnum(r11,s11)+1);
        M(nnum(r11,s11),nnum(r11,s11)+m)    = -MP(nnum(r11,s11),nnum(r11,s11)+m);
        M(nnum(r11,s11),nnum(r11,s11)-m)    = -MP(nnum(r11,s11),nnum(r11,s11)-m);
    end
    r11 = p1+p;
    nnum(r11,s11)=(r11-1)*m+s11;
    MP(nnum(r11,s11),nnum(r11,s11)-1)   = +1/Rv1;
    MP(nnum(r11,s11),nnum(r11,s11)+1)   =  1/Rv1;
    MP(nnum(r11,s11),nnum(r11,s11)-m)   =  1/Rw1;
    M(nnum(r11,s11),nnum(r11,s11))      = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)-m);
    M(nnum(r11,s11),nnum(r11,s11)+1)    = -MP(nnum(r11,s11),nnum(r11,s11)+1);
    M(nnum(r11,s11),nnum(r11,s11)-1)    = -MP(nnum(r11,s11),nnum(r11,s11)-1);
    M(nnum(r11,s11),nnum(r11,s11)-m)    = -MP(nnum(r11,s11),nnum(r11,s11)-m);
    end
end

%% stator teeth:1,3,5,7,9,11

% Columns left sides of teeth r11 >= p1
for dent = 1:6
    s11 = mal+1+(dent-1)*(2*mds+mes+ma);
    r11= p1;
    nnum(r11,s11)=(r11-1)*m+s11;
    MP(nnum(r11,s11),nnum(r11,s11)-1) = 1/Rvag;
    MP(nnum(r11,s11),nnum(r11,s11)+1) = 1/(Rvag);
    MP(nnum(r11,s11),nnum(r11,s11)-m) = 1/(Rwag);
    MP(nnum(r11,s11),nnum(r11,s11)+m) = 1/(Rw20/2+Rwag/2);
    M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
    M(nnum(r11,s11),nnum(r11,s11)-1)    = -MP(nnum(r11,s11),nnum(r11,s11)-1);
    M(nnum(r11,s11),nnum(r11,s11)+1)    = -MP(nnum(r11,s11),nnum(r11,s11)+1);
    M(nnum(r11,s11),nnum(r11,s11)+m)    = -MP(nnum(r11,s11),nnum(r11,s11)+m);
    M(nnum(r11,s11),nnum(r11,s11)-m)    = -MP(nnum(r11,s11),nnum(r11,s11)-m);
    r11 = p1+1;
    nnum(r11,s11)=(r11-1)*m+s11;
    MP(nnum(r11,s11),nnum(r11,s11)+1) = 1/Rv20;
    MP(nnum(r11,s11),nnum(r11,s11)-1) = 1/(Rv20/2+Rv1/2);
    MP(nnum(r11,s11),nnum(r11,s11)-m) = 1/(Rw20/2+Rwag/2);
    MP(nnum(r11,s11),nnum(r11,s11)+m) = 1/Rw20;
    M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
    M(nnum(r11,s11),nnum(r11,s11)-1)    = -MP(nnum(r11,s11),nnum(r11,s11)-1);
    M(nnum(r11,s11),nnum(r11,s11)+1)    = -MP(nnum(r11,s11),nnum(r11,s11)+1);
    M(nnum(r11,s11),nnum(r11,s11)+m)    = -MP(nnum(r11,s11),nnum(r11,s11)+m);
    M(nnum(r11,s11),nnum(r11,s11)-m)    = -MP(nnum(r11,s11),nnum(r11,s11)-m);
    for r11 = p1+2:p1+p-1
        nnum(r11,s11)=(r11-1)*m+s11;
        MP(nnum(r11,s11),nnum(r11,s11)+1) = 1/Rv20;
        MP(nnum(r11,s11),nnum(r11,s11)-1) = 1/(Rv20/2+Rv1/2);
        MP(nnum(r11,s11),nnum(r11,s11)-m) = 1/Rw20;
        MP(nnum(r11,s11),nnum(r11,s11)+m) = 1/Rw20;
        M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
        M(nnum(r11,s11),nnum(r11,s11)-1)    = -MP(nnum(r11,s11),nnum(r11,s11)-1);
        M(nnum(r11,s11),nnum(r11,s11)+1)    = -MP(nnum(r11,s11),nnum(r11,s11)+1);
        M(nnum(r11,s11),nnum(r11,s11)+m)    = -MP(nnum(r11,s11),nnum(r11,s11)+m);
        M(nnum(r11,s11),nnum(r11,s11)-m)    = -MP(nnum(r11,s11),nnum(r11,s11)-m);
    end
    r11 = p1+p;
    nnum(r11,s11)=(r11-1)*m+s11;
    MP(nnum(r11,s11),nnum(r11,s11)-1)   = +1/(Rv20/2+Rv1/2);
    MP(nnum(r11,s11),nnum(r11,s11)+1)   =  1/Rv20;
    MP(nnum(r11,s11),nnum(r11,s11)-m)   =  1/Rw20;
    M(nnum(r11,s11),nnum(r11,s11))      = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)-m);
    M(nnum(r11,s11),nnum(r11,s11)+1)    = -MP(nnum(r11,s11),nnum(r11,s11)+1);
    M(nnum(r11,s11),nnum(r11,s11)-1)    = -MP(nnum(r11,s11),nnum(r11,s11)-1);
    M(nnum(r11,s11),nnum(r11,s11)-m)    = -MP(nnum(r11,s11),nnum(r11,s11)-m);
end

% Columns right sides r11 >= p1
for dent = 1:6
    s11 = mal+mds+(dent-1)*(2*mds+mes+ma);
    r11= p1;
    nnum(r11,s11)=(r11-1)*m+s11;
    MP(nnum(r11,s11),nnum(r11,s11)-1) = 1/Rvag;
    MP(nnum(r11,s11),nnum(r11,s11)+1) = 1/(Rvag);
    MP(nnum(r11,s11),nnum(r11,s11)-m) = 1/(Rwag);
    MP(nnum(r11,s11),nnum(r11,s11)+m) = 1/(Rw20/2+Rwag/2);
    M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
    M(nnum(r11,s11),nnum(r11,s11)-1)    = -MP(nnum(r11,s11),nnum(r11,s11)-1);
    M(nnum(r11,s11),nnum(r11,s11)+1)    = -MP(nnum(r11,s11),nnum(r11,s11)+1);
    M(nnum(r11,s11),nnum(r11,s11)+m)    = -MP(nnum(r11,s11),nnum(r11,s11)+m);
    M(nnum(r11,s11),nnum(r11,s11)-m)    = -MP(nnum(r11,s11),nnum(r11,s11)-m);
    r11 = p1+1;
    nnum(r11,s11)=(r11-1)*m+s11;
    MP(nnum(r11,s11),nnum(r11,s11)-1) = 1/Rv20;
    MP(nnum(r11,s11),nnum(r11,s11)+1) = 1/(Rv20/2+Rv1/2);
    MP(nnum(r11,s11),nnum(r11,s11)-m) = 1/(Rw20/2+Rwag/2);
    MP(nnum(r11,s11),nnum(r11,s11)+m) = 1/Rw20;
    M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
    M(nnum(r11,s11),nnum(r11,s11)-1)    = -MP(nnum(r11,s11),nnum(r11,s11)-1);
    M(nnum(r11,s11),nnum(r11,s11)+1)    = -MP(nnum(r11,s11),nnum(r11,s11)+1);
    M(nnum(r11,s11),nnum(r11,s11)+m)    = -MP(nnum(r11,s11),nnum(r11,s11)+m);
    M(nnum(r11,s11),nnum(r11,s11)-m)    = -MP(nnum(r11,s11),nnum(r11,s11)-m);
    for r11 = p1+2:p1+10
        nnum(r11,s11)=(r11-1)*m+s11;
        MP(nnum(r11,s11),nnum(r11,s11)-1) = 1/Rv20;
        MP(nnum(r11,s11),nnum(r11,s11)+1) = 1/(Rv20/2+Rv1/2);
        MP(nnum(r11,s11),nnum(r11,s11)-m) = 1/Rw20;
        MP(nnum(r11,s11),nnum(r11,s11)+m) = 1/Rw20;
        M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
        M(nnum(r11,s11),nnum(r11,s11)-1)    = -MP(nnum(r11,s11),nnum(r11,s11)-1);
        M(nnum(r11,s11),nnum(r11,s11)+1)    = -MP(nnum(r11,s11),nnum(r11,s11)+1);
        M(nnum(r11,s11),nnum(r11,s11)+m)    = -MP(nnum(r11,s11),nnum(r11,s11)+m);
        M(nnum(r11,s11),nnum(r11,s11)-m)    = -MP(nnum(r11,s11),nnum(r11,s11)-m);
    end
    for r11 = p1+11:p1+p-1
        nnum(r11,s11)=(r11-1)*m+s11;
        MP(nnum(r11,s11),nnum(r11,s11)-1) = 1/Rv20;
        MP(nnum(r11,s11),nnum(r11,s11)+1) = 1/(Rv20);
        MP(nnum(r11,s11),nnum(r11,s11)-m) = 1/Rw20;
        MP(nnum(r11,s11),nnum(r11,s11)+m) = 1/Rw20;
        M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
        M(nnum(r11,s11),nnum(r11,s11)-1)    = -MP(nnum(r11,s11),nnum(r11,s11)-1);
        M(nnum(r11,s11),nnum(r11,s11)+1)    = -MP(nnum(r11,s11),nnum(r11,s11)+1);
        M(nnum(r11,s11),nnum(r11,s11)+m)    = -MP(nnum(r11,s11),nnum(r11,s11)+m);
        M(nnum(r11,s11),nnum(r11,s11)-m)    = -MP(nnum(r11,s11),nnum(r11,s11)-m);
    end
    r11 = p1+p;
    nnum(r11,s11)=(r11-1)*m+s11;
    MP(nnum(r11,s11),nnum(r11,s11)-1)   =  1/(Rv20);
    MP(nnum(r11,s11),nnum(r11,s11)+1)   =  1/Rv20;
    MP(nnum(r11,s11),nnum(r11,s11)-m)   =  1/Rw20;
    M(nnum(r11,s11),nnum(r11,s11))      = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)-m);
    M(nnum(r11,s11),nnum(r11,s11)+1)    = -MP(nnum(r11,s11),nnum(r11,s11)+1);
    M(nnum(r11,s11),nnum(r11,s11)-1)    = -MP(nnum(r11,s11),nnum(r11,s11)-1);
    M(nnum(r11,s11),nnum(r11,s11)-m)    = -MP(nnum(r11,s11),nnum(r11,s11)-m);
end

% the inside nodes of these teeth
for dent = 1:6
    for s11 = mal+2+(dent-1)*(2*mds+mes+ma):mal+mds-1+(dent-1)*(2*mds+mes+ma)
    r11= p1;
    nnum(r11,s11)=(r11-1)*m+s11;
    MP(nnum(r11,s11),nnum(r11,s11)-1) = 1/Rvag;
    MP(nnum(r11,s11),nnum(r11,s11)+1) = 1/Rvag;
    MP(nnum(r11,s11),nnum(r11,s11)-m) = 1/(Rwag);
    MP(nnum(r11,s11),nnum(r11,s11)+m) = 1/(Rw20/2+Rwag/2);
    M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
    M(nnum(r11,s11),nnum(r11,s11)-1)    = -MP(nnum(r11,s11),nnum(r11,s11)-1);
    M(nnum(r11,s11),nnum(r11,s11)+1)    = -MP(nnum(r11,s11),nnum(r11,s11)+1);
    M(nnum(r11,s11),nnum(r11,s11)+m)    = -MP(nnum(r11,s11),nnum(r11,s11)+m);
    M(nnum(r11,s11),nnum(r11,s11)-m)    = -MP(nnum(r11,s11),nnum(r11,s11)-m);
    r11 = p1+1;
    nnum(r11,s11)=(r11-1)*m+s11;
    MP(nnum(r11,s11),nnum(r11,s11)-1) = 1/Rv20;
    MP(nnum(r11,s11),nnum(r11,s11)+1) = 1/(Rv20);
    MP(nnum(r11,s11),nnum(r11,s11)-m) = 1/(Rw20/2+Rwag/2);
    MP(nnum(r11,s11),nnum(r11,s11)+m) = 1/Rw20;
    M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
    M(nnum(r11,s11),nnum(r11,s11)-1)    = -MP(nnum(r11,s11),nnum(r11,s11)-1);
    M(nnum(r11,s11),nnum(r11,s11)+1)    = -MP(nnum(r11,s11),nnum(r11,s11)+1);
    M(nnum(r11,s11),nnum(r11,s11)+m)    = -MP(nnum(r11,s11),nnum(r11,s11)+m);
    M(nnum(r11,s11),nnum(r11,s11)-m)    = -MP(nnum(r11,s11),nnum(r11,s11)-m);
    for r11 = p1+2:p1+p-1
        nnum(r11,s11)=(r11-1)*m+s11;
        MP(nnum(r11,s11),nnum(r11,s11)-1) = 1/Rv20;
        MP(nnum(r11,s11),nnum(r11,s11)+1) = 1/(Rv20);
        MP(nnum(r11,s11),nnum(r11,s11)-m) = 1/Rw20;
        MP(nnum(r11,s11),nnum(r11,s11)+m) = 1/Rw20;
        M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
        M(nnum(r11,s11),nnum(r11,s11)-1)    = -MP(nnum(r11,s11),nnum(r11,s11)-1);
        M(nnum(r11,s11),nnum(r11,s11)+1)    = -MP(nnum(r11,s11),nnum(r11,s11)+1);
        M(nnum(r11,s11),nnum(r11,s11)+m)    = -MP(nnum(r11,s11),nnum(r11,s11)+m);
        M(nnum(r11,s11),nnum(r11,s11)-m)    = -MP(nnum(r11,s11),nnum(r11,s11)-m);
    end
    r11 = p1+p;
    nnum(r11,s11)=(r11-1)*m+s11;
    MP(nnum(r11,s11),nnum(r11,s11)-1)   = +1/Rv20;
    MP(nnum(r11,s11),nnum(r11,s11)+1)   =  1/Rv20;
    MP(nnum(r11,s11),nnum(r11,s11)-m)   =  1/Rw20;
    M(nnum(r11,s11),nnum(r11,s11))      = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)-m);
    M(nnum(r11,s11),nnum(r11,s11)+1)    = -MP(nnum(r11,s11),nnum(r11,s11)+1);
    M(nnum(r11,s11),nnum(r11,s11)-1)    = -MP(nnum(r11,s11),nnum(r11,s11)-1);
    M(nnum(r11,s11),nnum(r11,s11)-m)    = -MP(nnum(r11,s11),nnum(r11,s11)-m);
    end
end

%% stator teeth: 2,4,6,8,10,12

% Columns right sides of the teeth r11 >= p1
for dent = 1:6
    s11 = mal+(2*mds+mes)+(dent-1)*(2*mds+mes+ma);
    r11= p1;
    nnum(r11,s11)=(r11-1)*m+s11;
    MP(nnum(r11,s11),nnum(r11,s11)-1) = 1/Rvag;
    MP(nnum(r11,s11),nnum(r11,s11)+1) = 1/(Rvag);
    MP(nnum(r11,s11),nnum(r11,s11)-m) = 1/(Rwag);
    MP(nnum(r11,s11),nnum(r11,s11)+m) = 1/(Rw20/2+Rwag/2);
    M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
    M(nnum(r11,s11),nnum(r11,s11)-1)    = -MP(nnum(r11,s11),nnum(r11,s11)-1);
    M(nnum(r11,s11),nnum(r11,s11)+1)    = -MP(nnum(r11,s11),nnum(r11,s11)+1);
    M(nnum(r11,s11),nnum(r11,s11)+m)    = -MP(nnum(r11,s11),nnum(r11,s11)+m);
    M(nnum(r11,s11),nnum(r11,s11)-m)    = -MP(nnum(r11,s11),nnum(r11,s11)-m);
    r11 = p1+1;
    nnum(r11,s11)=(r11-1)*m+s11;
    MP(nnum(r11,s11),nnum(r11,s11)-1) = 1/Rv20;
    MP(nnum(r11,s11),nnum(r11,s11)+1) = 1/(Rv20/2+Rv1/2);
    MP(nnum(r11,s11),nnum(r11,s11)-m) = 1/(Rw20/2+Rwag/2);
    MP(nnum(r11,s11),nnum(r11,s11)+m) = 1/Rw20;
    M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
    M(nnum(r11,s11),nnum(r11,s11)-1)    = -MP(nnum(r11,s11),nnum(r11,s11)-1);
    M(nnum(r11,s11),nnum(r11,s11)+1)    = -MP(nnum(r11,s11),nnum(r11,s11)+1);
    M(nnum(r11,s11),nnum(r11,s11)+m)    = -MP(nnum(r11,s11),nnum(r11,s11)+m);
    M(nnum(r11,s11),nnum(r11,s11)-m)    = -MP(nnum(r11,s11),nnum(r11,s11)-m);
    for r11 = p1+2:p1+p-1
        nnum(r11,s11)=(r11-1)*m+s11;
        MP(nnum(r11,s11),nnum(r11,s11)-1) = 1/Rv20;
        MP(nnum(r11,s11),nnum(r11,s11)+1) = 1/(Rv20/2+Rv1/2);
        MP(nnum(r11,s11),nnum(r11,s11)-m) = 1/Rw20;
        MP(nnum(r11,s11),nnum(r11,s11)+m) = 1/Rw20;
        M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
        M(nnum(r11,s11),nnum(r11,s11)-1)    = -MP(nnum(r11,s11),nnum(r11,s11)-1);
        M(nnum(r11,s11),nnum(r11,s11)+1)    = -MP(nnum(r11,s11),nnum(r11,s11)+1);
        M(nnum(r11,s11),nnum(r11,s11)+m)    = -MP(nnum(r11,s11),nnum(r11,s11)+m);
        M(nnum(r11,s11),nnum(r11,s11)-m)    = -MP(nnum(r11,s11),nnum(r11,s11)-m);
    end
    r11 = p1+p;
    nnum(r11,s11)=(r11-1)*m+s11;
    MP(nnum(r11,s11),nnum(r11,s11)+1)   = +1/(Rv20/2+Rv1/2);
    MP(nnum(r11,s11),nnum(r11,s11)-1)   =  1/Rv20;
    MP(nnum(r11,s11),nnum(r11,s11)-m)   =  1/Rw20;
    M(nnum(r11,s11),nnum(r11,s11))      = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)-m);
    M(nnum(r11,s11),nnum(r11,s11)+1)    = -MP(nnum(r11,s11),nnum(r11,s11)+1);
    M(nnum(r11,s11),nnum(r11,s11)-1)    = -MP(nnum(r11,s11),nnum(r11,s11)-1);
    M(nnum(r11,s11),nnum(r11,s11)-m)    = -MP(nnum(r11,s11),nnum(r11,s11)-m);
end

% Columns left sides of the teeth r11 >= p1
for dent = 1:6
    s11 = mal+mds+mes+1+(dent-1)*(2*mds+mes+ma);
    r11= p1;
    nnum(r11,s11)=(r11-1)*m+s11;
    MP(nnum(r11,s11),nnum(r11,s11)-1) = 1/Rvag;
    MP(nnum(r11,s11),nnum(r11,s11)+1) = 1/(Rvag);
    MP(nnum(r11,s11),nnum(r11,s11)-m) = 1/(Rwag);
    MP(nnum(r11,s11),nnum(r11,s11)+m) = 1/(Rw20/2+Rwag/2);
    M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
    M(nnum(r11,s11),nnum(r11,s11)-1)    = -MP(nnum(r11,s11),nnum(r11,s11)-1);
    M(nnum(r11,s11),nnum(r11,s11)+1)    = -MP(nnum(r11,s11),nnum(r11,s11)+1);
    M(nnum(r11,s11),nnum(r11,s11)+m)    = -MP(nnum(r11,s11),nnum(r11,s11)+m);
    M(nnum(r11,s11),nnum(r11,s11)-m)    = -MP(nnum(r11,s11),nnum(r11,s11)-m);
    r11 = p1+1;
    nnum(r11,s11)=(r11-1)*m+s11;
    MP(nnum(r11,s11),nnum(r11,s11)+1) = 1/Rv20;
    MP(nnum(r11,s11),nnum(r11,s11)-1) = 1/(Rv20/2+Rv1/2);
    MP(nnum(r11,s11),nnum(r11,s11)-m) = 1/(Rw20/2+Rwag/2);
    MP(nnum(r11,s11),nnum(r11,s11)+m) = 1/Rw20;
    M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
    M(nnum(r11,s11),nnum(r11,s11)-1)    = -MP(nnum(r11,s11),nnum(r11,s11)-1);
    M(nnum(r11,s11),nnum(r11,s11)+1)    = -MP(nnum(r11,s11),nnum(r11,s11)+1);
    M(nnum(r11,s11),nnum(r11,s11)+m)    = -MP(nnum(r11,s11),nnum(r11,s11)+m);
    M(nnum(r11,s11),nnum(r11,s11)-m)    = -MP(nnum(r11,s11),nnum(r11,s11)-m);
    for r11 = p1+2:p1+10
        nnum(r11,s11)=(r11-1)*m+s11;
        MP(nnum(r11,s11),nnum(r11,s11)+1) = 1/Rv20;
        MP(nnum(r11,s11),nnum(r11,s11)-1) = 1/(Rv20/2+Rv1/2);
        MP(nnum(r11,s11),nnum(r11,s11)-m) = 1/Rw20;
        MP(nnum(r11,s11),nnum(r11,s11)+m) = 1/Rw20;
        M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
        M(nnum(r11,s11),nnum(r11,s11)-1)    = -MP(nnum(r11,s11),nnum(r11,s11)-1);
        M(nnum(r11,s11),nnum(r11,s11)+1)    = -MP(nnum(r11,s11),nnum(r11,s11)+1);
        M(nnum(r11,s11),nnum(r11,s11)+m)    = -MP(nnum(r11,s11),nnum(r11,s11)+m);
        M(nnum(r11,s11),nnum(r11,s11)-m)    = -MP(nnum(r11,s11),nnum(r11,s11)-m);
    end
    for r11 = p1+11:p1+p-1
        nnum(r11,s11)=(r11-1)*m+s11;
        MP(nnum(r11,s11),nnum(r11,s11)-1) = 1/Rv20;
        MP(nnum(r11,s11),nnum(r11,s11)+1) = 1/(Rv20);
        MP(nnum(r11,s11),nnum(r11,s11)-m) = 1/Rw20;
        MP(nnum(r11,s11),nnum(r11,s11)+m) = 1/Rw20;
        M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
        M(nnum(r11,s11),nnum(r11,s11)-1)    = -MP(nnum(r11,s11),nnum(r11,s11)-1);
        M(nnum(r11,s11),nnum(r11,s11)+1)    = -MP(nnum(r11,s11),nnum(r11,s11)+1);
        M(nnum(r11,s11),nnum(r11,s11)+m)    = -MP(nnum(r11,s11),nnum(r11,s11)+m);
        M(nnum(r11,s11),nnum(r11,s11)-m)    = -MP(nnum(r11,s11),nnum(r11,s11)-m);
    end
    r11 = p1+p;
    nnum(r11,s11)=(r11-1)*m+s11;
    MP(nnum(r11,s11),nnum(r11,s11)-1)   =  1/(Rv20);
    MP(nnum(r11,s11),nnum(r11,s11)+1)   =  1/Rv20;
    MP(nnum(r11,s11),nnum(r11,s11)-m)   =  1/Rw20;
    M(nnum(r11,s11),nnum(r11,s11))      = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)-m);
    M(nnum(r11,s11),nnum(r11,s11)+1)    = -MP(nnum(r11,s11),nnum(r11,s11)+1);
    M(nnum(r11,s11),nnum(r11,s11)-1)    = -MP(nnum(r11,s11),nnum(r11,s11)-1);
    M(nnum(r11,s11),nnum(r11,s11)-m)    = -MP(nnum(r11,s11),nnum(r11,s11)-m);
end

% the inside nodes of these teeth
for dent = 1:6
    for s11 = mal+mes+mds+2+(dent-1)*(2*mds+mes+ma):mal+mes+2*mds-1+(dent-1)*(2*mds+mes+ma)
    r11= p1;
    nnum(r11,s11)=(r11-1)*m+s11;
    MP(nnum(r11,s11),nnum(r11,s11)-1) = 1/Rvag;
    MP(nnum(r11,s11),nnum(r11,s11)+1) = 1/(Rvag);
    MP(nnum(r11,s11),nnum(r11,s11)-m) = 1/(Rwag);
    MP(nnum(r11,s11),nnum(r11,s11)+m) = 1/(Rw20/2+Rwag/2);
    M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
    M(nnum(r11,s11),nnum(r11,s11)-1)    = -MP(nnum(r11,s11),nnum(r11,s11)-1);
    M(nnum(r11,s11),nnum(r11,s11)+1)    = -MP(nnum(r11,s11),nnum(r11,s11)+1);
    M(nnum(r11,s11),nnum(r11,s11)+m)    = -MP(nnum(r11,s11),nnum(r11,s11)+m);
    M(nnum(r11,s11),nnum(r11,s11)-m)    = -MP(nnum(r11,s11),nnum(r11,s11)-m);
    r11 = p1+1;
    nnum(r11,s11)=(r11-1)*m+s11;
    MP(nnum(r11,s11),nnum(r11,s11)-1) = 1/Rv20;
    MP(nnum(r11,s11),nnum(r11,s11)+1) = 1/(Rv20);
    MP(nnum(r11,s11),nnum(r11,s11)-m) = 1/(Rw20/2+Rwag/2);
    MP(nnum(r11,s11),nnum(r11,s11)+m) = 1/Rw20;
    M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
    M(nnum(r11,s11),nnum(r11,s11)-1)    = -MP(nnum(r11,s11),nnum(r11,s11)-1);
    M(nnum(r11,s11),nnum(r11,s11)+1)    = -MP(nnum(r11,s11),nnum(r11,s11)+1);
    M(nnum(r11,s11),nnum(r11,s11)+m)    = -MP(nnum(r11,s11),nnum(r11,s11)+m);
    M(nnum(r11,s11),nnum(r11,s11)-m)    = -MP(nnum(r11,s11),nnum(r11,s11)-m);
    for r11 = p1+2:p1+p-1
        nnum(r11,s11)=(r11-1)*m+s11;
        MP(nnum(r11,s11),nnum(r11,s11)-1) = 1/Rv20;
        MP(nnum(r11,s11),nnum(r11,s11)+1) = 1/(Rv20);
        MP(nnum(r11,s11),nnum(r11,s11)-m) = 1/Rw20;
        MP(nnum(r11,s11),nnum(r11,s11)+m) = 1/Rw20;
        M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
        M(nnum(r11,s11),nnum(r11,s11)-1)    = -MP(nnum(r11,s11),nnum(r11,s11)-1);
        M(nnum(r11,s11),nnum(r11,s11)+1)    = -MP(nnum(r11,s11),nnum(r11,s11)+1);
        M(nnum(r11,s11),nnum(r11,s11)+m)    = -MP(nnum(r11,s11),nnum(r11,s11)+m);
        M(nnum(r11,s11),nnum(r11,s11)-m)    = -MP(nnum(r11,s11),nnum(r11,s11)-m);
    end
    r11 = p1+p;
    nnum(r11,s11)=(r11-1)*m+s11;
    MP(nnum(r11,s11),nnum(r11,s11)-1)   = +1/Rv20;
    MP(nnum(r11,s11),nnum(r11,s11)+1)   =  1/Rv20;
    MP(nnum(r11,s11),nnum(r11,s11)-m)   =  1/Rw20;
    M(nnum(r11,s11),nnum(r11,s11))      = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)-m);
    M(nnum(r11,s11),nnum(r11,s11)+1)    = -MP(nnum(r11,s11),nnum(r11,s11)+1);
    M(nnum(r11,s11),nnum(r11,s11)-1)    = -MP(nnum(r11,s11),nnum(r11,s11)-1);
    M(nnum(r11,s11),nnum(r11,s11)-m)    = -MP(nnum(r11,s11),nnum(r11,s11)-m);
    end
end


%% stator slots
% left borders of the slots
for enc = 1:6
    s11 = mal+mds+1+(enc-1)*(2*mds+mes+ma);
    r11 = p1;
    nnum(r11,s11)=(r11-1)*m+s11;
    MP(nnum(r11,s11),nnum(r11,s11)-1) = 1/Rvag;
    MP(nnum(r11,s11),nnum(r11,s11)+1) = 1/Rvag;
    MP(nnum(r11,s11),nnum(r11,s11)-m) = 1/Rwag;
    MP(nnum(r11,s11),nnum(r11,s11)+m) = 1/(Rw1/2+Rwag/2);
    M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
    M(nnum(r11,s11),nnum(r11,s11)-1)  = -MP(nnum(r11,s11),nnum(r11,s11)-1);
    M(nnum(r11,s11),nnum(r11,s11)+1)  = -MP(nnum(r11,s11),nnum(r11,s11)+1);
    M(nnum(r11,s11),nnum(r11,s11)+m)  = -MP(nnum(r11,s11),nnum(r11,s11)+m);
    M(nnum(r11,s11),nnum(r11,s11)-m)  = -MP(nnum(r11,s11),nnum(r11,s11)-m);
    r11 = p1+1;
    nnum(r11,s11)=(r11-1)*m+s11;
    MP(nnum(r11,s11),nnum(r11,s11)-1) = 1/(Rv1/2+Rv20/2);
    MP(nnum(r11,s11),nnum(r11,s11)+1) = 1/Rv1;
    MP(nnum(r11,s11),nnum(r11,s11)+m) = 1/Rw1;
    MP(nnum(r11,s11),nnum(r11,s11)-m) = 1/(Rw1/2+Rwag/2);
    M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
    M(nnum(r11,s11),nnum(r11,s11)-1)  = -MP(nnum(r11,s11),nnum(r11,s11)-1);
    M(nnum(r11,s11),nnum(r11,s11)+1)  = -MP(nnum(r11,s11),nnum(r11,s11)+1);
    M(nnum(r11,s11),nnum(r11,s11)+m)  = -MP(nnum(r11,s11),nnum(r11,s11)+m);
    M(nnum(r11,s11),nnum(r11,s11)-m)  = -MP(nnum(r11,s11),nnum(r11,s11)-m);
    for r11 = p1+2:p1+9
        nnum(r11,s11)=(r11-1)*m+s11;
        MP(nnum(r11,s11),nnum(r11,s11)-1) = 1/(Rv20/2+Rv1/2);
        MP(nnum(r11,s11),nnum(r11,s11)+1) = 1/Rv1;
        MP(nnum(r11,s11),nnum(r11,s11)-m) = 1/Rw1;
        MP(nnum(r11,s11),nnum(r11,s11)+m) = 1/Rw1;
        M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
        M(nnum(r11,s11),nnum(r11,s11)-1)  = -MP(nnum(r11,s11),nnum(r11,s11)-1);
        M(nnum(r11,s11),nnum(r11,s11)+1)  = -MP(nnum(r11,s11),nnum(r11,s11)+1);
        M(nnum(r11,s11),nnum(r11,s11)+m)  = -MP(nnum(r11,s11),nnum(r11,s11)+m);
        M(nnum(r11,s11),nnum(r11,s11)-m)  = -MP(nnum(r11,s11),nnum(r11,s11)-m);
    end
    r11 = p1+10;
    nnum(r11,s11)=(r11-1)*m+s11;
    MP(nnum(r11,s11),nnum(r11,s11)-1) = 1/(Rv20/2+Rv1/2);
    MP(nnum(r11,s11),nnum(r11,s11)+1) = 1/Rv1;
    MP(nnum(r11,s11),nnum(r11,s11)-m) = 1/Rw1;
    MP(nnum(r11,s11),nnum(r11,s11)+m) = 1/(Rw1/2+Rw20/2);
    M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
    M(nnum(r11,s11),nnum(r11,s11)-1)  = -MP(nnum(r11,s11),nnum(r11,s11)-1);
    M(nnum(r11,s11),nnum(r11,s11)+1)  = -MP(nnum(r11,s11),nnum(r11,s11)+1);
    M(nnum(r11,s11),nnum(r11,s11)+m)  = -MP(nnum(r11,s11),nnum(r11,s11)+m);
    M(nnum(r11,s11),nnum(r11,s11)-m)  = -MP(nnum(r11,s11),nnum(r11,s11)-m);
    r11 = p1+11;
    nnum(r11,s11)=(r11-1)*m+s11;
    MP(nnum(r11,s11),nnum(r11,s11)-1) = 1/(Rv20);
    MP(nnum(r11,s11),nnum(r11,s11)+1) = 1/Rv20;
    MP(nnum(r11,s11),nnum(r11,s11)+m) = 1/Rw20;
    MP(nnum(r11,s11),nnum(r11,s11)-m) = 1/(Rw1/2+Rw20/2);
    M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
    M(nnum(r11,s11),nnum(r11,s11)-1)  = -MP(nnum(r11,s11),nnum(r11,s11)-1);
    M(nnum(r11,s11),nnum(r11,s11)+1)  = -MP(nnum(r11,s11),nnum(r11,s11)+1);
    M(nnum(r11,s11),nnum(r11,s11)+m)  = -MP(nnum(r11,s11),nnum(r11,s11)+m);
    M(nnum(r11,s11),nnum(r11,s11)-m)  = -MP(nnum(r11,s11),nnum(r11,s11)-m);
    for r11 = p1+12:p1+p-1
        nnum(r11,s11)=(r11-1)*m+s11;
        MP(nnum(r11,s11),nnum(r11,s11)-1) = 1/(Rv20);
        MP(nnum(r11,s11),nnum(r11,s11)+1) = 1/Rv20;
        MP(nnum(r11,s11),nnum(r11,s11)-m) = 1/Rw20;
        MP(nnum(r11,s11),nnum(r11,s11)+m) = 1/Rw20;
        M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
        M(nnum(r11,s11),nnum(r11,s11)-1)  = -MP(nnum(r11,s11),nnum(r11,s11)-1);
        M(nnum(r11,s11),nnum(r11,s11)+1)  = -MP(nnum(r11,s11),nnum(r11,s11)+1);
        M(nnum(r11,s11),nnum(r11,s11)+m)  = -MP(nnum(r11,s11),nnum(r11,s11)+m);
        M(nnum(r11,s11),nnum(r11,s11)-m)  = -MP(nnum(r11,s11),nnum(r11,s11)-m);
    end
    r11 = p1+p;
    nnum(r11,s11)=(r11-1)*m+s11;
    MP(nnum(r11,s11),nnum(r11,s11)-1)   = +1/Rv20;
    MP(nnum(r11,s11),nnum(r11,s11)+1)   =  1/Rv20;
    MP(nnum(r11,s11),nnum(r11,s11)-m)   =  1/Rw20;
    M(nnum(r11,s11),nnum(r11,s11))      = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)-m);
    M(nnum(r11,s11),nnum(r11,s11)+1)    = -MP(nnum(r11,s11),nnum(r11,s11)+1);
    M(nnum(r11,s11),nnum(r11,s11)-1)    = -MP(nnum(r11,s11),nnum(r11,s11)-1);
    M(nnum(r11,s11),nnum(r11,s11)-m)    = -MP(nnum(r11,s11),nnum(r11,s11)-m);
end

% right borders of the slots
for enc = 1:6
    s11 = mal+mds+mes+(enc-1)*(2*mds+mes+ma);
    r11 = p1;
    nnum(r11,s11)=(r11-1)*m+s11;
    MP(nnum(r11,s11),nnum(r11,s11)-1) = 1/Rvag;
    MP(nnum(r11,s11),nnum(r11,s11)+1) = 1/Rvag;
    MP(nnum(r11,s11),nnum(r11,s11)-m) = 1/Rwag;
    MP(nnum(r11,s11),nnum(r11,s11)+m) = 1/(Rw1/2+Rwag/2);
    M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
    M(nnum(r11,s11),nnum(r11,s11)-1)  = -MP(nnum(r11,s11),nnum(r11,s11)-1);
    M(nnum(r11,s11),nnum(r11,s11)+1)  = -MP(nnum(r11,s11),nnum(r11,s11)+1);
    M(nnum(r11,s11),nnum(r11,s11)+m)  = -MP(nnum(r11,s11),nnum(r11,s11)+m);
    M(nnum(r11,s11),nnum(r11,s11)-m)  = -MP(nnum(r11,s11),nnum(r11,s11)-m);
    r11 = p1+1;
    nnum(r11,s11)=(r11-1)*m+s11;
    MP(nnum(r11,s11),nnum(r11,s11)-1) = 1/Rv1;
    MP(nnum(r11,s11),nnum(r11,s11)+1) = 1/(Rv1/2+Rv20/2);
    MP(nnum(r11,s11),nnum(r11,s11)+m) = 1/Rw1;
    MP(nnum(r11,s11),nnum(r11,s11)-m) = 1/(Rw1/2+Rwag/2);
    M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
    M(nnum(r11,s11),nnum(r11,s11)-1)  = -MP(nnum(r11,s11),nnum(r11,s11)-1);
    M(nnum(r11,s11),nnum(r11,s11)+1)  = -MP(nnum(r11,s11),nnum(r11,s11)+1);
    M(nnum(r11,s11),nnum(r11,s11)+m)  = -MP(nnum(r11,s11),nnum(r11,s11)+m);
    M(nnum(r11,s11),nnum(r11,s11)-m)  = -MP(nnum(r11,s11),nnum(r11,s11)-m);
    for r11 = p1+2:p1+9
        nnum(r11,s11)=(r11-1)*m+s11;
        MP(nnum(r11,s11),nnum(r11,s11)+1) = 1/(Rv20/2+Rv1/2);
        MP(nnum(r11,s11),nnum(r11,s11)-1) = 1/Rv1;
        MP(nnum(r11,s11),nnum(r11,s11)-m) = 1/Rw1;
        MP(nnum(r11,s11),nnum(r11,s11)+m) = 1/Rw1;
        M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
        M(nnum(r11,s11),nnum(r11,s11)-1)  = -MP(nnum(r11,s11),nnum(r11,s11)-1);
        M(nnum(r11,s11),nnum(r11,s11)+1)  = -MP(nnum(r11,s11),nnum(r11,s11)+1);
        M(nnum(r11,s11),nnum(r11,s11)+m)  = -MP(nnum(r11,s11),nnum(r11,s11)+m);
        M(nnum(r11,s11),nnum(r11,s11)-m)  = -MP(nnum(r11,s11),nnum(r11,s11)-m);
    end
    r11 = p1+10;
    nnum(r11,s11)=(r11-1)*m+s11;
    MP(nnum(r11,s11),nnum(r11,s11)+1) = 1/(Rv20/2+Rv1/2);
    MP(nnum(r11,s11),nnum(r11,s11)-1) = 1/Rv1;
    MP(nnum(r11,s11),nnum(r11,s11)-m) = 1/Rw1;
    MP(nnum(r11,s11),nnum(r11,s11)+m) = 1/(Rw1/2+Rw20/2);
    M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
    M(nnum(r11,s11),nnum(r11,s11)-1)  = -MP(nnum(r11,s11),nnum(r11,s11)-1);
    M(nnum(r11,s11),nnum(r11,s11)+1)  = -MP(nnum(r11,s11),nnum(r11,s11)+1);
    M(nnum(r11,s11),nnum(r11,s11)+m)  = -MP(nnum(r11,s11),nnum(r11,s11)+m);
    M(nnum(r11,s11),nnum(r11,s11)-m)  = -MP(nnum(r11,s11),nnum(r11,s11)-m);
    r11 = p1+11;
    nnum(r11,s11)=(r11-1)*m+s11;
    MP(nnum(r11,s11),nnum(r11,s11)-1) = 1/(Rv20);
    MP(nnum(r11,s11),nnum(r11,s11)+1) = 1/Rv20;
    MP(nnum(r11,s11),nnum(r11,s11)+m) = 1/Rw20;
    MP(nnum(r11,s11),nnum(r11,s11)-m) = 1/(Rw1/2+Rw20/2);
    M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
    M(nnum(r11,s11),nnum(r11,s11)-1)  = -MP(nnum(r11,s11),nnum(r11,s11)-1);
    M(nnum(r11,s11),nnum(r11,s11)+1)  = -MP(nnum(r11,s11),nnum(r11,s11)+1);
    M(nnum(r11,s11),nnum(r11,s11)+m)  = -MP(nnum(r11,s11),nnum(r11,s11)+m);
    M(nnum(r11,s11),nnum(r11,s11)-m)  = -MP(nnum(r11,s11),nnum(r11,s11)-m);
    for r11 = p1+12:p1+p-1
        nnum(r11,s11)=(r11-1)*m+s11;
        MP(nnum(r11,s11),nnum(r11,s11)-1) = 1/(Rv20);
        MP(nnum(r11,s11),nnum(r11,s11)+1) = 1/Rv20;
        MP(nnum(r11,s11),nnum(r11,s11)-m) = 1/Rw20;
        MP(nnum(r11,s11),nnum(r11,s11)+m) = 1/Rw20;
        M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
        M(nnum(r11,s11),nnum(r11,s11)-1)  = -MP(nnum(r11,s11),nnum(r11,s11)-1);
        M(nnum(r11,s11),nnum(r11,s11)+1)  = -MP(nnum(r11,s11),nnum(r11,s11)+1);
        M(nnum(r11,s11),nnum(r11,s11)+m)  = -MP(nnum(r11,s11),nnum(r11,s11)+m);
        M(nnum(r11,s11),nnum(r11,s11)-m)  = -MP(nnum(r11,s11),nnum(r11,s11)-m);
    end
    r11 = p1+p;
    nnum(r11,s11)=(r11-1)*m+s11;
    MP(nnum(r11,s11),nnum(r11,s11)-1)   = +1/Rv20;
    MP(nnum(r11,s11),nnum(r11,s11)+1)   =  1/Rv20;
    MP(nnum(r11,s11),nnum(r11,s11)-m)   =  1/Rw20;
    M(nnum(r11,s11),nnum(r11,s11))      = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)-m);
    M(nnum(r11,s11),nnum(r11,s11)+1)    = -MP(nnum(r11,s11),nnum(r11,s11)+1);
    M(nnum(r11,s11),nnum(r11,s11)-1)    = -MP(nnum(r11,s11),nnum(r11,s11)-1);
    M(nnum(r11,s11),nnum(r11,s11)-m)    = -MP(nnum(r11,s11),nnum(r11,s11)-m);
end


% inner nodes of the rotor slots
for enc = 1:6
    for s11 = mal+mds+2+(enc-1)*(2*mds+mes+ma):mal+mds+mes-1+(enc-1)*(2*mds+mes+ma)
        r11 = p1;
        nnum(r11,s11)=(r11-1)*m+s11;
        MP(nnum(r11,s11),nnum(r11,s11)-1) = 1/Rvag;
        MP(nnum(r11,s11),nnum(r11,s11)+1) = 1/Rvag;
        MP(nnum(r11,s11),nnum(r11,s11)-m) = 1/Rwag;
        MP(nnum(r11,s11),nnum(r11,s11)+m) = 1/(Rw1/2+Rwag/2);
        M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
        M(nnum(r11,s11),nnum(r11,s11)-1)  = -MP(nnum(r11,s11),nnum(r11,s11)-1);
        M(nnum(r11,s11),nnum(r11,s11)+1)  = -MP(nnum(r11,s11),nnum(r11,s11)+1);
        M(nnum(r11,s11),nnum(r11,s11)+m)  = -MP(nnum(r11,s11),nnum(r11,s11)+m);
        M(nnum(r11,s11),nnum(r11,s11)-m)  = -MP(nnum(r11,s11),nnum(r11,s11)-m);
        r11 = p1+1;
        nnum(r11,s11)=(r11-1)*m+s11;
        MP(nnum(r11,s11),nnum(r11,s11)-1) = 1/Rv1;
        MP(nnum(r11,s11),nnum(r11,s11)+1) = 1/(Rv1);
        MP(nnum(r11,s11),nnum(r11,s11)+m) = 1/Rw1;
        MP(nnum(r11,s11),nnum(r11,s11)-m) = 1/(Rw1/2+Rwag/2);
        M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
        M(nnum(r11,s11),nnum(r11,s11)-1)  = -MP(nnum(r11,s11),nnum(r11,s11)-1);
        M(nnum(r11,s11),nnum(r11,s11)+1)  = -MP(nnum(r11,s11),nnum(r11,s11)+1);
        M(nnum(r11,s11),nnum(r11,s11)+m)  = -MP(nnum(r11,s11),nnum(r11,s11)+m);
        M(nnum(r11,s11),nnum(r11,s11)-m)  = -MP(nnum(r11,s11),nnum(r11,s11)-m);
        for r11 = p1+2:p1+9
            nnum(r11,s11)=(r11-1)*m+s11;
            MP(nnum(r11,s11),nnum(r11,s11)+1) = 1/(Rv1);
            MP(nnum(r11,s11),nnum(r11,s11)-1) = 1/Rv1;
            MP(nnum(r11,s11),nnum(r11,s11)-m) = 1/Rw1;
            MP(nnum(r11,s11),nnum(r11,s11)+m) = 1/Rw1;
            M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
            M(nnum(r11,s11),nnum(r11,s11)-1)  = -MP(nnum(r11,s11),nnum(r11,s11)-1);
            M(nnum(r11,s11),nnum(r11,s11)+1)  = -MP(nnum(r11,s11),nnum(r11,s11)+1);
            M(nnum(r11,s11),nnum(r11,s11)+m)  = -MP(nnum(r11,s11),nnum(r11,s11)+m);
            M(nnum(r11,s11),nnum(r11,s11)-m)  = -MP(nnum(r11,s11),nnum(r11,s11)-m);
        end
        r11 = p1+10;
        nnum(r11,s11)=(r11-1)*m+s11;
        MP(nnum(r11,s11),nnum(r11,s11)+1) = 1/(Rv1);
        MP(nnum(r11,s11),nnum(r11,s11)-1) = 1/Rv1;
        MP(nnum(r11,s11),nnum(r11,s11)-m) = 1/Rw1;
        MP(nnum(r11,s11),nnum(r11,s11)+m) = 1/(Rw1/2+Rw20/2);
        M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
        M(nnum(r11,s11),nnum(r11,s11)-1)  = -MP(nnum(r11,s11),nnum(r11,s11)-1);
        M(nnum(r11,s11),nnum(r11,s11)+1)  = -MP(nnum(r11,s11),nnum(r11,s11)+1);
        M(nnum(r11,s11),nnum(r11,s11)+m)  = -MP(nnum(r11,s11),nnum(r11,s11)+m);
        M(nnum(r11,s11),nnum(r11,s11)-m)  = -MP(nnum(r11,s11),nnum(r11,s11)-m);
        r11 = p1+11;
        nnum(r11,s11)=(r11-1)*m+s11;
        MP(nnum(r11,s11),nnum(r11,s11)-1) = 1/(Rv20);
        MP(nnum(r11,s11),nnum(r11,s11)+1) = 1/Rv20;
        MP(nnum(r11,s11),nnum(r11,s11)+m) = 1/Rw20;
        MP(nnum(r11,s11),nnum(r11,s11)-m) = 1/(Rw1/2+Rw20/2);
        M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
        M(nnum(r11,s11),nnum(r11,s11)-1)  = -MP(nnum(r11,s11),nnum(r11,s11)-1);
        M(nnum(r11,s11),nnum(r11,s11)+1)  = -MP(nnum(r11,s11),nnum(r11,s11)+1);
        M(nnum(r11,s11),nnum(r11,s11)+m)  = -MP(nnum(r11,s11),nnum(r11,s11)+m);
        M(nnum(r11,s11),nnum(r11,s11)-m)  = -MP(nnum(r11,s11),nnum(r11,s11)-m);
        for r11 = p1+12:p1+p-1
            nnum(r11,s11)=(r11-1)*m+s11;
            MP(nnum(r11,s11),nnum(r11,s11)-1) = 1/(Rv20);
            MP(nnum(r11,s11),nnum(r11,s11)+1) = 1/Rv20;
            MP(nnum(r11,s11),nnum(r11,s11)-m) = 1/Rw20;
            MP(nnum(r11,s11),nnum(r11,s11)+m) = 1/Rw20;
            M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
            M(nnum(r11,s11),nnum(r11,s11)-1)  = -MP(nnum(r11,s11),nnum(r11,s11)-1);
            M(nnum(r11,s11),nnum(r11,s11)+1)  = -MP(nnum(r11,s11),nnum(r11,s11)+1);
            M(nnum(r11,s11),nnum(r11,s11)+m)  = -MP(nnum(r11,s11),nnum(r11,s11)+m);
            M(nnum(r11,s11),nnum(r11,s11)-m)  = -MP(nnum(r11,s11),nnum(r11,s11)-m);
        end
        r11 = p1+p;
        nnum(r11,s11)=(r11-1)*m+s11;
        MP(nnum(r11,s11),nnum(r11,s11)-1)   = +1/Rv20;
        MP(nnum(r11,s11),nnum(r11,s11)+1)   =  1/Rv20;
        MP(nnum(r11,s11),nnum(r11,s11)-m)   =  1/Rw20;
        M(nnum(r11,s11),nnum(r11,s11))      = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)-m);
        M(nnum(r11,s11),nnum(r11,s11)+1)    = -MP(nnum(r11,s11),nnum(r11,s11)+1);
        M(nnum(r11,s11),nnum(r11,s11)-1)    = -MP(nnum(r11,s11),nnum(r11,s11)-1);
        M(nnum(r11,s11),nnum(r11,s11)-m)    = -MP(nnum(r11,s11),nnum(r11,s11)-m);
    end
end

%% excitation vector E
E = zeros(1,nn);
%nnum2 = zeros(30,353);

for dent = 1:6
    s11 = mal+(dent-1)*(2*mds+mes+ma);
    for r11 = p1+1:p1+p
        nnum2(r11,s11) = (r11-1)*m+s11;
        E(nnum2(r11,s11)) = (-1)^dent*Br0*la*h/2;
        E(nnum2(r11,s11)+1) = E(nnum2(r11,s11));
    end
end

for dent = 1:6
    s11 = mal+2*mds+mes+(dent-1)*(2*mds+mes+ma);
    for r11 = p1+1:p1+p
        nnum2(r11,s11) = (r11-1)*m+s11;
        E(nnum2(r11,s11)) = (-1)^dent*Br0*la*h/2;
        E(nnum2(r11,s11)+1) = E(nnum2(r11,s11));
    end
end
E(nn) = 0;

%% solve the system
C=M\E';
%% Post processing
%% Flux through phases
Ups = zeros(p+p1,m);
nnum3 = zeros(p+p1,m);

for i=1:(p+p1)
    for j=1:m
        nnum3(i,j)=(i-1)*m+j;
        Ups(i,j) = C(nnum3(i,j));        
    end
end

rangeTooth1 = mal+1:mal+mds;
rangeTooth2 = mal+mds+mes+1:mal+2*mds+mes;
rangeTooth3 = mal+(2*mds+mes+ma)+1:mal+(2*mds+mes+ma)+mds;
rangeTooth4 = mal+(2*mds+mes+ma)+mds+mes+1:mal+(2*mds+mes+ma)+2*mds+mes;
rangeTooth5 = mal+2*(2*mds+mes+ma)+1:mal+2*(2*mds+mes+ma)+mds;
rangeTooth6 = mal+2*(2*mds+mes+ma)+mes+mds+1:mal+2*(2*mds+mes+ma)+mes+2*mds;
rangeTooth7 = mal+3*(2*mds+mes+ma)+1:mal+3*(2*mds+mes+ma)+mds;
rangeTooth8 = mal+3*(2*mds+mes+ma)+mes+mds+1:mal+3*(2*mds+mes+ma)+mes+2*mds;
rangeTooth9 = mal+4*(2*mds+mes+ma)+1:mal+4*(2*mds+mes+ma)+mds;
rangeTooth10 = mal+4*(2*mds+mes+ma)+mds+mes+1:mal+4*(2*mds+mes+ma)+2*mds+mes;
rangeTooth11 = mal+5*(2*mds+mes+ma)+1:mal+5*(2*mds+mes+ma)+mds;
rangeTooth12 = mal+5*(2*mds+mes+ma)+mds+mes+1:mal+5*(2*mds+mes+ma)+2*mds+mes;

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
%% Magnetic Induction
% Calculation of magnetic flux density in the air gap
for s11 = 2:m-1
    nnumB1   = (p1-1)*m+s11;
    Bx(s11,1) = mu0*(C(nnumB1-1)-C(nnumB1+1))/2/l;
end

nnumB1   = (p1-1)*m+1;
Bx(1,1) = mu0*(-C(nnumB1-1+m)-C(nnumB1+1))/2/l;
Bx(m,1) = mu0*(C(nnumB1-1+m-1)+C(nnumB1))/2/l;

for i = 1:m
    nnumB31   = (p1)*m+i;
    nnumB22   = (p1-1)*m+i;
    nnumB23   = (p1-2)*m+i;
    nnumB24   = (p1-3)*m+i;
    By1(i,1) = (C(nnumB23)-C(nnumB22))*MP(nnumB23,nnumB22)/(l*la);
    By2(i,1) = (C(nnumB22)-C(nnumB31))*MP(nnumB22,nnumB31)/(l*la);
    By3(i,1) = (C(nnumB24)-C(nnumB23))*MP(nnumB24,nnumB23)/(l*la);
    By(i,1)  = (By1(i,1)+By3(i,1)+By2(i,1))/3;
end

FluxA = zeros(ni,1);
FluxB = zeros(ni,1);
FluxC = zeros(ni,1);

FluxB(1)     = sum(Phi(:,2))+sum(Phi(:,3))+sum(Phi(:,8))+sum(Phi(:,9));
FluxA(1)     = sum(Phi(:,1))+sum(Phi(:,12))+sum(Phi(:,6))+sum(Phi(:,7));
FluxC(1)     = sum(Phi(:,4))+sum(Phi(:,5))+sum(Phi(:,10))+sum(Phi(:,11));

toc

move_rotor

save FluxA_720_sym_e_15_MeanRadius FluxA;
save FluxB_720_sym_e_15_MeanRadius FluxB;
save FluxC_720_sym_e_15_MeanRadius FluxC;
save B_720_sym_e_15_MeanRadius Bx By;
save Vars_720_sym_e_15_MeanRadius l u m la R Rint Rext Rm;
toc