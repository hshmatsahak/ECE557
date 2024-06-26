function xdot = cartpend(in1,u)
%CARTPEND
%    XDOT = CARTPEND(IN1,U)

%    This function was generated by the Symbolic Math Toolbox version 9.3.
%    05-Dec-2023 19:14:35

theta = in1(3,:);
theta_dot = in1(4,:);
z_dot = in1(2,:);
t2 = cos(theta);
t3 = sin(theta);
t4 = theta_dot.^2;
t6 = z_dot.*(2.1e+1./2.0);
t8 = u.*1.7265;
t5 = t3.^2;
mt1 = [z_dot;-(t6-t8-t2.*t3.*2.254+t3.*t4.*7.5946e-2)./(t5.*(2.3e+1./1.0e+2)+1.1911);theta_dot];
mt2 = [-(t3.*(-1.392678e+1)+t2.*(t6-t8)+t2.*t3.*t4.*7.5946e-2)./(t5.*7.5946e-2+3.9330122e-1)];
xdot = [mt1;mt2];
end
