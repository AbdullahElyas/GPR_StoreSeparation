clc;
clear;

% 6dof simulation

% load the model data and GPR model
GPR = RegressionModel(0.1,0.9); %split train test dataset
regression_model = GPR.GPR_Model;

% 6dof data for Mach 0.95 AOA 0
data = GPR.SimData(3).SixDOF_Processed_Data;

% Store Parameters
IXX= 27.116;
IYY  = 488.094;
IZZ  = 488.094;
MOI_Matrix = [IXX 0 0; 0 IYY 0; 0 0 IZZ];
mass_store = 907.185; %store mass

AOA_aircraft = 0; % AOA in degrees
Beta_aircraft = 0; % Beta in degrees
Mach_aircraft = 0.95; % Mach number
flowAltitude = 26000/3.28; %meters
[~,sound_speed,~,~] = atmosisa(flowAltitude);
V_airspeed = Mach_aircraft*sound_speed;

% wind to body axis transformation matrix
wind2body = [cos(AOA_aircraft) 0 -sin(AOA_aircraft);...
    0 1 0;...
    sin(AOA_aircraft) 0 cos(AOA_aircraft)]*...
    [cos(Beta_aircraft) sin(Beta_aircraft) 0;...
    -sin(Beta_aircraft) cos(Beta_aircraft) 0;...
    0 0 1];
V_body = wind2body*[V_airspeed 0 0]';

% Gavitational force
g = 9.81; % m/s^2
F_Gravity = [0 0 mass_store*g]';


% define the initial conditions
% create a zero initial state vector
x0 = zeros(1,14);
x0(1) = 0; % initial bodyx position
x0(2) = 0; % initial bodyy position
x0(3) = 0; % initial bodyz position
x0(4) = 0; % initial roll angle in radian
x0(5) = 0; % initial pitch angle in radian
x0(6) = 0; % initial yaw angle in radian
x0(7) = V_body(1); % initial x velocity
x0(8) = V_body(2); % initial y velocity
x0(9) = V_body(3); % initial z velocity
x0(10) = 0; % initial roll rate
x0(11) = 0; % initial pitch rate
x0(12) = 0; % initial yaw rate
x0(13) = 0; % initial AOA store in radian
x0(14) = 0; % initial Beta store in radian
x0(15) = 0; % initial inertialx position
x0(16) = 0; % initial inertialy position
x0(17) = 0; % initial inertialz position

% define the time vector
tspan = [0:0.0005:0.5];

% solve the ODE using ode45(RK4 Method with adaptive time step)
[t,x] = ode45(@(t,x) rigid_body_dynamics(t,x,F_Gravity,mass_store,MOI_Matrix,regression_model,V_airspeed,Mach_aircraft,AOA_aircraft),tspan,x0);

% plot
figure
plot(t,x(:,15))
xlabel('Time (s)')
ylabel('X position')
hold on
plot(data.time,-data.X)
legend('Simulation','Data')
figure
plot(t,x(:,16))
xlabel('Time (s)')
ylabel('Y position')
hold on
plot(data.time,-data.Z)
legend('Simulation','Data')
figure
plot(t,x(:,17))
xlabel('Time (s)')

ylabel('Z position')
hold on
plot(data.time,-data.Y)

legend('Simulation','Data')
figure
plot(t,x(:,4)*57.3)
xlabel('Time (s)')
ylabel('Roll angle (deg)')
hold on
plot(data.time,data.Euler_Roll)
legend('Simulation','Data')

figure
plot(t,x(:,5)*57.3)
xlabel('Time (s)')
ylabel('Pitch angle (deg)')
hold on
plot(data.time,data.Euler_Pitch)


legend('Simulation','Data')
figure
plot(t,x(:,6)*57.3)
xlabel('Time (s)')
ylabel('Yaw angle (deg)')

hold on
plot(data.time,data.Euler_Yaw)
legend('Simulation','Data')




% create a function representing the rigid body dynamics
function [dx] = rigid_body_dynamics(t,x,Force_gravity,mass_store,MOI_Matrix,regression_model,V_airspeed,Mach_aircraft,AOA_aircraft)
% define the state variables
% x = [x y z phi theta psi u v w p q r]
% define the constants
m = mass_store; % mass
I = MOI_Matrix; % inertia matrix
% define the state variables
xb = x(1);
yb = x(2);
zb = x(3);
phi = x(4);
theta = x(5);
psi = x(6);
u = x(7);
v = x(8);
w = x(9);
p = x(10);
q = x(11);
r = x(12);

% define the rotation matrix from body to inertial axis
R = [cos(theta)*cos(psi) cos(theta)*sin(psi) -sin(theta);...
    sin(phi)*sin(theta)*cos(psi)-cos(phi)*sin(psi) sin(phi)*sin(theta)*sin(psi)+cos(phi)*cos(psi) sin(phi)*cos(theta);...
    cos(phi)*sin(theta)*cos(psi)+sin(phi)*sin(psi) cos(phi)*sin(theta)*sin(psi)-sin(phi)*cos(psi) cos(phi)*cos(theta)];

% determine axis of rotation and angle of rotation from rotation matrix
t
rotationVector = rotmat2vec3d(R);
thetaZ_sim = -rotationVector(2);
thetaY_sim = -rotationVector(3);
thetaX_sim = -rotationVector(1);
% estimate the aero body forces and moments
[force_body_aero_sim,moment_body_aero_sim] = GPR_predict(x,Mach_aircraft,AOA_aircraft,regression_model);

% sim to body axis covention xbody = -xsim. ybody = -zsim, zbody = -ysim
force_body_aero = [-force_body_aero_sim(1) -force_body_aero_sim(3) -force_body_aero_sim(2)];
moment_body_aero = [-moment_body_aero_sim(1) -moment_body_aero_sim(3) -moment_body_aero_sim(2)];
% convert the force_gravity to body axis
force_gravity = R*Force_gravity;

% ejector force and moment
[total_ejector_force_inertial,total_ejector_moment_inertial] = ejector_force(x,thetaZ_sim);
total_ejector_force_body = R*total_ejector_force_inertial';
total_ejector_moment_body = R*total_ejector_moment_inertial';

% define the total body forces and moments

% F_body = F_aero + F_gravity + F_ejector
force_body = force_body_aero + force_gravity'+total_ejector_force_body';
% M_body = M_aero + M_ejector
moment_body = moment_body_aero+total_ejector_moment_body';

% define the linear and angular accelerations
dx(1) = u;
dx(2) = v;
dx(3) = w;
dx(4) = p + q*sin(phi)*tan(theta) + r*cos(phi)*tan(theta);
dx(5) = q*cos(phi) - r*sin(phi);
dx(6) = q*sin(phi)/cos(theta) + r*cos(phi)/cos(theta);
dx(7) = force_body(1)/m - q*w + r*v;
dx(8) = force_body(2)/m - r*u + p*w;
dx(9) = force_body(3)/m - p*v + q*u;
dx(10) = (I(2,2)-I(3,3))*q*r/I(1,1) + moment_body(1)/I(1,1);
dx(11) = (I(3,3)-I(1,1))*p*r/I(2,2) + moment_body(2)/I(2,2);
dx(12) = (I(1,1)-I(2,2))*p*q/I(3,3) + moment_body(3)/I(3,3);

% extimate AOA and Beta
alpha = atan(w/u);
beta = asin(v/sqrt(u^2+v^2+w^2));
% define the state vector
dx(13) = alpha;
dx(14) = beta;
% inertial velocities xe, ye, ze
dx(15) = u*cos(theta)*cos(psi) + v*(sin(phi)*sin(theta)*cos(psi)-cos(phi)*sin(psi)) + w*(cos(phi)*sin(theta)*cos(psi)+sin(phi)*sin(psi))-V_airspeed;
dx(16) = u*cos(theta)*sin(psi) + v*(sin(phi)*sin(theta)*sin(psi)+cos(phi)*cos(psi)) + w*(cos(phi)*sin(theta)*sin(psi)-sin(phi)*cos(psi));
dx(17) = -u*sin(theta) + v*sin(phi)*cos(theta) + w*cos(phi)*cos(theta);
dx = dx';
end

% create a function to estimate the body forces and moments
function [Force_body_aero,Moment_body_aero] = GPR_predict(x,Mach_aircraft,AOA_aircraft,regression_model)
INPUT = zeros(1,5);
INPUT(1) = -x(14); % x position converting back to sim axis
INPUT(2) = -x(16); % y position converting back to sim axis
INPUT(3) = -x(15); % z position converting back to sim axis
INPUT(4) = x(13)*57.3; % AOA store in deg
INPUT(5) = x(14)*57.3; % Beta store in deg
INPUT(6) = Mach_aircraft;
INPUT(7) = AOA_aircraft;
% define the forces and moments
force_body_x_pred =  predict(regression_model.Mdl_force_body_x,INPUT);
force_body_y_pred =  predict(regression_model.Mdl_force_body_y,INPUT);
force_body_z_pred =  predict(regression_model.Mdl_force_body_z,INPUT);
moment_body_x_pred =  predict(regression_model.Mdl_moment_body_x,INPUT);
moment_body_y_pred =  predict(regression_model.Mdl_moment_body_y,INPUT);
moment_body_z_pred =  predict(regression_model.Mdl_moment_body_z,INPUT);
% define the body forces and moments
Force_body_aero = [force_body_x_pred force_body_y_pred force_body_z_pred];
Moment_body_aero = [moment_body_x_pred moment_body_y_pred moment_body_z_pred];

% assume damping derivatives for force and moments and add to the forces and moments
% Force_body_aero = Force_body_aero + [-0.1*x(7) -0.1*x(8) -0.1*x(9)];
Moment_body_aero = Moment_body_aero + [0.1*3000*x(10)*0 -0.1*5000*x(11)*0 -0.1*89000*x(12)*0];
end

% Function for Ejector Force and Moment (Stroke Length 0.100584m)
function [total_ejector_force_inertial,total_ejector_moment_inertial] = ejector_force(x,thetaZ_sim)
% define the ejector force
front_ejector_momentarm = 0.179832;
rear_ejector_momentarm = 0.329184;
theta = -thetaZ_sim;
dfront = abs(x(17) - front_ejector_momentarm*theta);
dback = abs(x(17) + rear_ejector_momentarm*theta);
% Initialize forces and moments
FY = 0;
MZ = 0;
% Check if front and rear ejector stroke length is less than 0.100584m
if dfront <= 0.100584
    FY = -10676.0;
    MZ = 1920.0;
end
if dback <= 0.100584
    FY = FY - 42703.0;
    MZ = MZ - 14057.0;
end
total_ejector_force_inertial = [0 0 -FY];
total_ejector_moment_inertial = [0 -MZ 0];
end










