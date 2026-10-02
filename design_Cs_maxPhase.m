function [Cs, a, b, info] = design_Cs_maxPhase(k_inf, omega_r, p, wpk, phi_deg, varargin)
%DESIGN_CS_MAXPHASE  Max-phase shaping filter under Eq_m constraints + phase cap.
%
%   C_s(s) = (sum_{k=0}^p a_k * (s/omega_r)^k) / (sum_{k=0}^p b_k * (s/omega_r)^k)
%   with a0=b0=1 and, for m=1..p:
%       Eq_m:  sum_{k+l=2m-2} a_k b_l (-1)^l  +  sum_{k+l=2m-1} a_k b_l (-1)^l = 0.
%
%   Inputs:
%     k_inf   = a_p / b_p  (high-frequency gain)
%     omega_r > 0          (rad/s)
%     p       >= 2         (order)
%     wpk     > 0          (rad/s), desired phase-peak frequency
%     phi_deg              target maximum phase at wpk (deg). Enforced as phase(wpk) <= phi_deg,
%                          and the optimizer tries to get as close as possible from below.
%
%   Output:
%     Cs    : tf object in s-domain (numerator/denominator in powers of s)
%     a, b  : coefficient vectors [a0 ... ap] and [b0 ... bp] (ASCENDING in k)
%     info  : diagnostics
%
%   Requires: Optimization Toolbox (fmincon).
%
%   Example:
%     [Cs,a,b,info] = design_Cs_maxPhase(0.3, 2*pi*10, 4, 2*pi*26, 40);

% -------------------- options --------------------
opt.NumStarts   = 30;      % multistart attempts
opt.EpsRel      = 1e-3;    % relative step around wpk for derivative approx
opt.GridDec     = 0.6;     % +/- decades around wpk for "peak dominance" constraints
opt.GridN       = 9;       % number of grid points in that band (including wpk)
opt.WnMin       = 1e-3;    % min dimensionless wn in x=s/omega_r domain
opt.ZetaMin     = 0.02;    % min damping (>0 ensures LHP)
opt.BetaMin     = 1e-3;    % min first-order time constant in x-domain
opt.LogBound    = 8;       % bounds on log-parameters
opt.Verbose     = true;
opt.PhiMarginDeg = 0.05;   % for example    

if ~isempty(varargin)
    opt = parse_name_value(opt, varargin{:});
end

% -------------------- input checks --------------------
validateattributes(k_inf,   {'double'},{'real','finite','scalar'});
validateattributes(omega_r, {'double'},{'real','finite','positive','scalar'});
validateattributes(p,       {'double'},{'real','finite','integer','scalar','>=',2});
validateattributes(wpk,     {'double'},{'real','finite','positive','scalar'});
validateattributes(phi_deg, {'double'},{'real','finite','scalar'});

phiT = deg2rad(phi_deg);

% factor counts (stable-by-construction parameterization)
n2 = floor(p/2);     % number of 2nd-order sections
n1 = mod(p,2);       % number of 1st-order sections (0 or 1)

% decision variables (unconstrained -> exp -> positive)
% ordering:
%   den: [log_beta_d (if n1) ; log_wn_d(1:n2) ; log_zeta_d(1:n2)]
%   num: [log_beta_n (if n1) ; log_wn_n(1:n2) ; log_zeta_n(1:n2)]
nd = n1 + 2*n2;
nx = 2*nd;

lb = -opt.LogBound * ones(nx,1);
ub =  opt.LogBound * ones(nx,1);

% fmincon options
fopts = optimoptions('fmincon',...
    'Algorithm','sqp',...
    'Display','none',...
    'MaxFunctionEvaluations',2e5,...
    'MaxIterations',2e3,...
    'ConstraintTolerance',1e-3,...
    'OptimalityTolerance',1e-3);

% best feasible / infeasible trackers
bestFeas = struct('z',[],'phi',-Inf,'a',[],'b',[],'exitflag',[],'output',[], ...
                  'c',[],'ceq',[],'viol',Inf,'fval',Inf);
bestInf  = struct('z',[],'phi',-Inf,'a',[],'b',[],'exitflag',[],'output',[], ...
                  'c',[],'ceq',[],'viol',Inf,'fval',Inf);

% -------------------- multistart loop --------------------
for trial = 1:opt.NumStarts
    z0 = initial_guess();

    problem = createOptimProblem('fmincon',...
        'x0',z0,'objective',@obj,'nonlcon',@nonlcon,'lb',lb,'ub',ub,'options',fopts);

    try
        [zsol,fval,exitflag,output] = fmincon(problem);
    catch
        break; % solver not available or error
    end

    [c,ceq,a_try,b_try,phi_try] = constraints_and_phase(zsol);
    viol = norm([max(c,0); ceq],2);
    isFeas = all(c <= 1e-6) && all(abs(ceq) <= 1e-6);

    if isFeas
        if phi_try > bestFeas.phi
            bestFeas = pack_best(zsol,phi_try,a_try,b_try,exitflag,output,c,ceq,viol,fval);
        end
    else
        % keep best "least-infeasible" as fallback
        if (viol < bestInf.viol) || (abs(viol-bestInf.viol) < 1e-9 && phi_try > bestInf.phi)
            bestInf = pack_best(zsol,phi_try,a_try,b_try,exitflag,output,c,ceq,viol,fval);
        end
    end
end

useFeas = ~isempty(bestFeas.z);
best = bestFeas;
if ~useFeas
    best = bestInf;
end

if isempty(best.z)
    error('design_Cs_maxPhase:NoSolution','No solution attempt completed (solver not available or failed immediately).');
end

a = best.a;  % ascending [a0..ap]
b = best.b;  % ascending [b0..bp]

% Build tf in s-domain:
%   sum a_k*(s/omega_r)^k = sum (a_k/omega_r^k) s^k
pow = omega_r .^ (0:p);
num_desc = fliplr(a ./ pow);   % descending in s
den_desc = fliplr(b ./ pow);

Cs = tf(num_desc, den_desc);

% info
info = struct();
info.feasible         = useFeas;
info.exitflag         = best.exitflag;
info.output           = best.output;
info.constraint_c     = best.c;
info.constraint_ceq   = best.ceq;
info.constraint_violation_L2 = best.viol;

info.bestPhase_rad    = best.phi;
info.bestPhase_deg    = rad2deg(best.phi);
info.targetPhase_deg  = phi_deg;
info.phaseError_deg   = info.bestPhase_deg - phi_deg;   % should be <=0 if feasible

info.k_inf_achieved   = a(end)/b(end);
info.a = a;
info.b = b;

if opt.Verbose
    if info.feasible
        fprintf('[design_Cs_maxPhase] FEASIBLE: p=%d, phase(wpk)=%.3f deg (<= %.3f), k_inf=%.6g\n',...
            p, info.bestPhase_deg, phi_deg, info.k_inf_achieved);
    else
        fprintf('[design_Cs_maxPhase] INFEASIBLE best: p=%d, phase(wpk)=%.3f deg, ||viol||2=%.3e, k_inf=%.6g\n',...
            p, info.bestPhase_deg, info.constraint_violation_L2, info.k_inf_achieved);
    end
end

% ==================== nested helpers ====================

    function z0 = initial_guess()
        % Heuristic: center natural frequencies around r = wpk/omega_r in x-domain.
        r = wpk/omega_r;

        wn_d = max(opt.WnMin, r * exp(0.7*randn(n2,1)));
        wn_n = max(opt.WnMin, r * exp(0.7*randn(n2,1)));

        zeta_d = max(opt.ZetaMin, 0.3 + 0.4*rand(n2,1));
        zeta_n = max(opt.ZetaMin, 0.3 + 0.4*rand(n2,1));

        if n1==1
            beta_d = max(opt.BetaMin, r * exp(0.7*randn(1,1)));
            beta_n = max(opt.BetaMin, r * exp(0.7*randn(1,1)));
            z0 = [log(beta_d - opt.BetaMin + 1e-12);
                  log(wn_d   - opt.WnMin   + 1e-12);
                  log(zeta_d - opt.ZetaMin + 1e-12);
                  log(beta_n - opt.BetaMin + 1e-12);
                  log(wn_n   - opt.WnMin   + 1e-12);
                  log(zeta_n - opt.ZetaMin + 1e-12)];
        else
            z0 = [log(wn_d   - opt.WnMin   + 1e-12);
                  log(zeta_d - opt.ZetaMin + 1e-12);
                  log(wn_n   - opt.WnMin   + 1e-12);
                  log(zeta_n - opt.ZetaMin + 1e-12)];
        end

        z0 = z0(:);
        if numel(z0) ~= nx
            error('Internal sizing mismatch (nx=%d, got %d).', nx, numel(z0));
        end
        z0 = min(max(z0,lb),ub);
    end

    function f = obj(z)
        % maximize phase at wpk (from below phiT due to constraint)
        [~,~,~,~,phi0] = constraints_and_phase(z);
        f = -phi0 + 1e-4*(z.'*z); % regularization
    end

    function [c,ceq] = nonlcon(z)
        [c,ceq] = constraints_and_phase(z);
    end

    function [c,ceq,a_here,b_here,phi0] = constraints_and_phase(z)
        [a_here,b_here] = build_ab_from_z(z);

        % ----- Eq_m constraints via convolution: coeffs of A(x)*B(-x)
        bt = b_here .* ((-1).^(0:p));      % b_l * (-1)^l
        conv_abt = conv(a_here, bt);       % n=0..2p, length 2p+1

        ceq_Eq = zeros(p,1);
        for m = 1:p
            nA = 2*m-2;
            nB = 2*m-1;
            ceq_Eq(m) = conv_abt(nA+1) + conv_abt(nB+1);
        end

        % ----- k_inf constraint
        ceq_k = a_here(end) - k_inf*b_here(end);

        % ----- phase peak derivative constraint at wpk
        epsRel = opt.EpsRel;
        w1 = wpk*(1-epsRel);
        w2 = wpk;
        w3 = wpk*(1+epsRel);

        phis = local_unwrap_phase([w1,w2,w3], a_here, b_here); % radians, unwrapped
        phis = shift_branch_to_target(phis, 2, phiT);
        phi_m = phis(1); phi0 = phis(2); phi_p = phis(3);

        dphi = (phi_p - phi_m) / (w3 - w1); % should be 0 at extremum

        % ----- "wpk is a local max" inequality constraints on a band
        wgrid = wpk * logspace(-opt.GridDec, opt.GridDec, opt.GridN);
        wgrid = unique([wgrid(:); wpk]);
        wgrid = sort(wgrid);

        phig = local_unwrap_phase(wgrid.', a_here, b_here);
        [~,i0] = min(abs(wgrid - wpk));
        phig = shift_branch_to_target(phig, i0, phiT);

        phi0g = phig(i0);
        phig_others = phig;
        phig_others(i0) = [];
        c_peakBand = phig_others - phi0g; % enforce <= 0

        % curvature negative (encourages a strict max)
        ddphi = (phi_p - 2*phi0 + phi_m) / ((w3-w2)^2);
        c_curv = ddphi; % <= 0

        % ----- maximum phase limit at wpk: phi0 <= phiT
        c_phiCap = phi0 - phiT; % <= 0

        c   = [c_peakBand(:); c_curv; c_phiCap];
        ceq = [ceq_Eq(:); ceq_k; dphi];
    end

    function [a_here,b_here] = build_ab_from_z(z)
        idx = 1;

        % denominator params
        if n1==1
            beta_d = opt.BetaMin + exp(z(idx)); idx = idx+1;
        else
            beta_d = [];
        end
        wn_d   = opt.WnMin   + exp(z(idx:idx+n2-1)); idx = idx+n2;
        zeta_d = opt.ZetaMin + exp(z(idx:idx+n2-1)); idx = idx+n2;

        % numerator params
        if n1==1
            beta_n = opt.BetaMin + exp(z(idx)); idx = idx+1;
        else
            beta_n = [];
        end
        wn_n   = opt.WnMin   + exp(z(idx:idx+n2-1)); idx = idx+n2;
        zeta_n = opt.ZetaMin + exp(z(idx:idx+n2-1)); idx = idx+n2;

        % build b(x) and a(x) in x = s/omega_r, ascending coefficients
        b_here = stable_poly_asc(p, beta_d, wn_d, zeta_d);
        a_here = stable_poly_asc(p, beta_n, wn_n, zeta_n);

        % enforce exact a0=b0=1
        b_here = b_here / b_here(1);
        a_here = a_here / a_here(1);
    end

    function polyAsc = stable_poly_asc(pLocal, beta1, wnVec, zetaVec)
        % polynomial in x with constant term 1, built from stable sections:
        % 1st order: (1 + x/beta)
        % 2nd order: (x/wn)^2 + 2*zeta*(x/wn) + 1
        polyAsc = 1; % ascending
        if ~isempty(beta1)
            polyAsc = conv(polyAsc, [1, 1/beta1]);
        end
        for i = 1:numel(wnVec)
            wn = wnVec(i);
            zt = zetaVec(i);
            polyAsc = conv(polyAsc, [1, 2*zt/wn, 1/(wn^2)]);
        end
        if numel(polyAsc) ~= (pLocal+1)
            error('Sectioning produced wrong degree (got %d, expected %d).', numel(polyAsc)-1, pLocal);
        end
    end

    function ph = local_unwrap_phase(wvec, aAsc, bAsc)
        % phase of C(jw) over wvec (assumed ordered), unwrapped
        Cjw = zeros(size(wvec));
        for ii = 1:numel(wvec)
            x = 1j*(wvec(ii)/omega_r);
            num = polyval(fliplr(aAsc), x); % fliplr => descending for polyval
            den = polyval(fliplr(bAsc), x);
            Cjw(ii) = num/den;
        end
        ph = unwrap(angle(Cjw));
    end

    function ph2 = shift_branch_to_target(ph, idx0, target)
        % Add an integer multiple of 2*pi so that ph(idx0) is closest to target.
        k = round((target - ph(idx0)) / (2*pi));
        ph2 = ph + 2*pi*k;
    end

    function s = pack_best(zsol,phi_try,a_try,b_try,exitflag,output,c,ceq,viol,fval)
        s = struct('z',zsol,'phi',phi_try,'a',a_try,'b',b_try,'exitflag',exitflag,...
                   'output',output,'c',c,'ceq',ceq,'viol',viol,'fval',fval);
    end

end

% -------- simple name-value parser --------
function opt = parse_name_value(opt, varargin)
if mod(numel(varargin),2)~=0
    error('Name-value pairs required.');
end
for i=1:2:numel(varargin)
    name = varargin{i};
    val  = varargin{i+1};
    if ~isfield(opt,name)
        error('Unknown option "%s".', name);
    end
    opt.(name) = val;
end
end