function [position, velocity, force, propeller, forceDerivative, t, costValue, elapsedTime, freeFalling, expectedForce] = solve_propelled_mass(inputStruct)

m = inputStruct.m;
g = inputStruct.g;

dT = inputStruct.dT;
T = inputStruct.T;

x0 = inputStruct.x0;
v0 = inputStruct.v0;
f0 = inputStruct.f0;

cost_multiplier = inputStruct.cost_multiplier;

M_fdot = inputStruct.M_fdot;
epsilon_relaxed = inputStruct.epsilon_relaxed;
K_dynamic = inputStruct.K_dynamic;
epsilon_dynamical = inputStruct.epsilon_dynamical;
K_hyperbolic = inputStruct.K_hyperbolic;
scaling_hyperbolic = inputStruct.scaling_hyperbolic;

%% Definition of the dynamics
X = casadi.MX.sym('X', 3);
U = casadi.MX.sym('U', 2);

x = X(1);
f = X(2);
v = X(3);
p = U(1);
f_dot = U(2);

X_dot = [v;
    f_dot;
    g + 1/m * (f + p)];

F = casadi.Function('continuous_dynamics', {X,U}, {X_dot});

%% Discretization
X_k1 = casadi.MX.sym('X_k1', 3);
X_k2 = casadi.MX.sym('X_k2', 3);
U_k1 = casadi.MX.sym('U_k1', 2);
U_k2 = casadi.MX.sym('U_k2', 2);

col_err = X_k1 + 0.5 * dT *  (F(X_k1, U_k1) + F(X_k2, U_k2)); %Implicit trapezoidal
%col_err = X_k1 + dT *  F(X_k1, U_k1); %explicit euler

F_k = casadi.Function('discrete_dynamics', {X_k1, U_k1, X_k2, U_k2}, {col_err});

%% Complementarity

rc = x * f - epsilon_relaxed;
dc = v * f + x * f_dot + K_dynamic * x * f - epsilon_dynamical;
delta = 1/cosh(scaling_hyperbolic * x);
hc = f_dot - delta * M_fdot + (1 - delta) * K_hyperbolic * f;

relaxed_complemenarity       = casadi.Function('Relaxed_Complemenarity', {X,U}, {rc});
dynamical_complemenarity  = casadi.Function('Dynamical_Complemenarity', {X,U}, {dc});
hyperbolic_complementarity  = casadi.Function('Hyperbolic_Complemenarity', {X,U}, {hc});

if strcmp(inputStruct.complementairity, 'Relaxed')
    complementarity = relaxed_complemenarity;
elseif strcmp(inputStruct.complementairity, 'Dynamical')
    complementarity = dynamical_complemenarity;
elseif strcmp(inputStruct.complementairity, 'Hyperbolic')
    complementarity = hyperbolic_complementarity;    
end

%% Problem solutions
N = round(T/dT);

t= 0 : dT : N*dT;

freeFalling = x0 + v0*t + 0.5*g*t.*t;
flying = freeFalling > 0;
freeFalling = flying .* freeFalling;
free_falling_velocity = v0 + g * t;
free_falling_velocity = flying .* free_falling_velocity;
in_contact = ~freeFalling;
expectedForce = m*abs(g)*in_contact;

opti = casadi.Opti();
Xo = opti.variable(3, N + 1);
Uo = opti.variable(2, N+1);

opti.subject_to(Xo(:,1) == [x0; f0; v0]); %initial state

for k = 1 : N
    opti.subject_to(Xo(:,k+1) == F_k(Xo(:,k),Uo(:,k), Xo(:,k+1),Uo(:,k+1))); %dynamics
    opti.subject_to(Xo(1,k) >= 0); %position positive
    opti.subject_to(Xo(2,k) >= 0); %force positive
    opti.subject_to(-M_fdot <= Uo(2, k)); %max force derivative
    opti.subject_to(Uo(2, k) <= M_fdot); %max force derivative
    opti.subject_to(complementarity(Xo(:,k),Uo(:,k)) <= 0) %complementarity
end

opti.subject_to(Xo(1,N+1) >= 0); %position positive
opti.subject_to(Xo(2,N+1) >= 0); %force positive
opti.subject_to(-M_fdot <= Uo(2, N+1)); %max force derivative
opti.subject_to(Uo(2, N+1) <= M_fdot); %max force derivative
opti.subject_to(complementarity(Xo(:,N+1),Uo(:,N+1)) <= 0) %complementarity
%opti.subject_to(Xo(2,N+1) == -m*g); %The force is equal to the weight at the end
%opti.subject_to(Uo(1,N) == 0); %No propeller at the end


opti.set_initial(Xo, [freeFalling; expectedForce; free_falling_velocity]);

cost = cost_multiplier*sumsqr(Uo(1,:));

opti.minimize(cost);


options = struct;
options.expand = true;
options.print_time = true;

options.ipopt.print_level = 0;
options.ipopt.linear_solver='ma57';
options.ipopt.ma57_pivtol = 1e-6;
options.ipopt.nlp_scaling_max_gradient = 100.0;
options.ipopt.nlp_scaling_min_value = 1e-6;
options.ipopt.tol = 1e-3;
options.ipopt.dual_inf_tol = 1000.0;
options.ipopt.compl_inf_tol = 1e-2;
options.ipopt.constr_viol_tol = 1e-4;
options.ipopt.acceptable_tol = 1e0;
options.ipopt.acceptable_iter =  2;
options.ipopt.acceptable_compl_inf_tol = 1.0;
options.ipopt.alpha_for_y = 'dual-and-full';
options.ipopt.max_iter =  4000;
options.ipopt.ma97_print_level = -10;
options.ipopt.warm_start_bound_frac = 1e-2;
options.ipopt.warm_start_bound_push = 1e-2;
options.ipopt.warm_start_mult_bound_push = 1e-2;
options.ipopt.warm_start_slack_bound_frac = 1e-2;
options.ipopt.warm_start_slack_bound_push = 1e-2;
options.ipopt.warm_start_init_point = 'yes';
options.ipopt.required_infeasibility_reduction = 0.8;
options.ipopt.perturb_dec_fact = 0.1;
options.ipopt.max_hessian_perturbation = 100.0;
options.ipopt.fast_step_computation = 'yes';
options.ipopt.hessian_approximation = 'limited-memory';

opti.solver('ipopt', options);

sol = opti.solve();

elapsedTime = sol.stats.t_wall_total;

xsol = sol.value(Xo);
usol = sol.value(Uo);
costValue = sol.value(cost)/cost_multiplier;

position = xsol(1,:);
force = xsol(2, :);
velocity = xsol(3, :);
propeller = usol(1, :);
forceDerivative = usol(2, :);

end