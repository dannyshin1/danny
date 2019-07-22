% Author :  Qida Yu                                                        
% This programe is implemented in matlab 2018a      
% Address:  Nanjing University of Aeronautics and Astronautics              

function [R_out, t_out] = OPnL(p1,p2,W1,W2)
%Thie version return all mimima to the end solution;
U=[W1,W2];
R = cat(3, rotx(pi/2), roty(pi/2), rotz(pi/2));
t = mean(U,2);

C_test = [];
t_test = [];

for i = 1:3
    % Make a random rotation
    pp1 = R(:,:,i) * (W1 - repmat(t, 1, size(W1,2)));
    pp2 = R(:,:,i) * (W2 - repmat(t, 1, size(W2,2)));
    %case1
    [C_est_i, t_est_i] = OPnL_main( p1,p2,pp1,pp2 );
    
    %rotate back to the original system
    for j = 1:size(t_est_i,2)
        t_est_i(:,j) = t_est_i(:,j) - C_est_i(:,:,j) * R(:,:,i) * t;
        C_est_i(:,:,j) = C_est_i(:,:,j) * R(:,:,i);
    end
    
	%stack all rotations and translations     
    if size(t_est_i,2) > 0
        C_test = cat(3,C_test, C_est_i);
        t_test = [t_test t_est_i];
    end   
end

R_out=C_test;
t_out=t_test;
end

function r = rotx(t)
% roty: rotation about x-axis
ct = cos(t);
st = sin(t);
r =    [1	0	0;
    0	ct	-st;
    0	st	ct];
end

function r = roty(t)
% roty: rotation about y-axis
ct = cos(t);
st = sin(t);
r =    [ct	0	st;
    0	1	0;
    -st	0	ct];
end

function r = rotz(t)
% rotz: rotation about z-axis
ct = cos(t);
st = sin(t);
r = [ct	-st	0
    st	ct	0
    0	0	1];

end
