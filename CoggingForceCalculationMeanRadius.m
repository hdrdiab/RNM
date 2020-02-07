load B_720_sym_e_15_MeanRadius;
load Vars_720_sym_e_15_MeanRadius

% Force calculation using Maxwell's tensor
mu0  = 4*pi*1e-7;        % air permeability (H/m)
for i=1:120
    for j=1:m
        ProdInd(i,j)  = Bx(j,i)*By(j,i);
    end
    torque(i)=0;
    for j=1:m-1
        % divided by 2 since we are averaging the ProdInd value over the
        % length between 2 nodes
        torque(i)  = torque(i)+(ProdInd(i,j+1)+ProdInd(i,j))*la*l/2/(mu0);
    end
end
Rm = (Rext+Rint)/2;
torque = torque*Rm;
figure
plot(1:length(torque),torque);

%save CoggForce_e15_MeanRadius torque;