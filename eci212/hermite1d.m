function hermite1d(L, numel, M0)

    if nargin < 3
        error('Usage: hermite1d(L, numel, M0)');
    end
    if L <= 0, error('L must be positive.'); end
    if numel < 1 || floor(numel) ~= numel, error('numel must be a positive integer.'); end

    if mod(numel,2) ~= 0
        fprintf('numel=%d is odd; increasing to %d so x=L/2 is a node.\n', numel, numel+1);
        numel = numel + 1;
    end

    EI = 1.0;              
    nn  = numel + 1;       
    ndof = 2*nn;           

    % Mesh
    x = linspace(0, L, nn).';
    le = L/numel;

    K = sparse(ndof, ndof);
    F = zeros(ndof, 1);

    % Assembly
    for e = 1:numel
        n1 = e;
        n2 = e+1;

        Ke = beamHermiteKe(EI, le);

        edofs = [2*n1-1, 2*n1, 2*n2-1, 2*n2];

        K(edofs, edofs) = K(edofs, edofs) + Ke;
    end

    midNode = numel/2 + 1;          
    thetaDof = 2*midNode;           
    F(thetaDof) = F(thetaDof) + M0; 

    fixed = [1, 2, 2*nn-1, 2*nn];      
    free  = setdiff(1:ndof, fixed);

    U = zeros(ndof, 1);
    U(free) = K(free, free) \ F(free);

    w     = U(1:2:end);
    theta = U(2:2:end);

    fprintf('\n--- hermite1d results ---\n');
    fprintf('L = %.6g, numel = %d, nn = %d, EI = %.6g\n', L, numel, nn, EI);
    fprintf('Applied moment M0 = %.6g at x = %.6g (node %d)\n', M0, x(midNode), midNode);
    fprintf('w(mid) = %.6e, theta(mid) = %.6e\n', w(midNode), theta(midNode));


    % ===== Plot FE vs exact =====
    % Exact solution uses the closed-form piecewise w(x) derived for a point moment at a=L/2.
    a = L/2;

    % Dense grid for exact curve
    xd = linspace(0, L, 1000).';
    wd = exact_w_clamped_clamped_point_moment(xd, L, EI, M0);

    % FE nodal curve (piecewise-linear between nodes for display)
    figure('Name','Hermite EB Beam: FE vs Exact (deflection)');
    plot(x, w, 'o-', 'LineWidth', 1.5); hold on;
    plot(xd, wd, '-', 'LineWidth', 2.0);
    grid on;
    xlabel('x'); ylabel('w(x)');
    title(sprintf('Clamped–clamped EB beam, point moment M_0 at x=L/2 (numel=%d)', numel));
    legend('FE (nodal)', 'Exact', 'Location', 'best');

    % Optional: rotations too (nodal vs exact)
    thetad = exact_theta_clamped_clamped_point_moment(xd, L, EI, M0);
    figure('Name','Hermite EB Beam: FE vs Exact (rotation)');
    plot(x, theta, 'o-', 'LineWidth', 1.5); hold on;
    plot(xd, thetad, '-', 'LineWidth', 2.0);
    grid on;
    xlabel('x'); ylabel('\theta(x)=w''(x)');
    title('Rotation: FE nodal vs exact');
    legend('FE (nodal)', 'Exact', 'Location', 'best');

% ===== add these local functions at the end of hermite1d.m (below beamHermiteKe) =====
function w = exact_w_clamped_clamped_point_moment(x, L, EI, M0)
%EXACT_W_CLAMPED_CLAMPED_POINT_MOMENT  Exact deflection for clamped-clamped beam
% with a concentrated moment M0 at x=L/2 (EI constant).
    a = L/2;
    w = zeros(size(x));

    left = (x <= a);
    xl = x(left);
    % 0 <= x <= L/2:
    w(left) = (M0/(8*EI))*xl.^2 - (M0/(4*EI*L))*xl.^3;

    right = ~left;
    xr = x(right);
    % L/2 <= x <= L:
    w(right) = (M0*L^2/(8*EI)) - (M0*L/(2*EI))*xr + (5*M0/(8*EI))*xr.^2 ...
               - (M0/(4*EI*L))*xr.^3;
end

function th = exact_theta_clamped_clamped_point_moment(x, L, EI, M0)
%EXACT_THETA_CLAMPED_CLAMPED_POINT_MOMENT  Exact rotation theta=w' for the same problem.
    a = L/2;
    th = zeros(size(x));

    left = (x <= a);
    xl = x(left);
    % 0 <= x <= L/2:
    th(left) = (M0/(4*EI))*xl - (3*M0/(4*EI*L))*xl.^2;

    right = ~left;
    xr = x(right);
    % L/2 <= x <= L:
    th(right) = -(M0*L/(2*EI)) + (5*M0/(4*EI))*xr - (3*M0/(4*EI*L))*xr.^2;
end









    % Plot deflection
%    figure('Name','Hermite EB Beam: deflection');
%    plot(x, w, 'o-', 'LineWidth', 1.5);
%    grid on;
%    xlabel('x');
%    ylabel('w(x)');
%    title(sprintf('Clamped–clamped EB beam, point moment M_0 at x=L/2 (numel=%d)', numel));
%
%    % Optional: show rotations too
%    figure('Name','Hermite EB Beam: rotation');
%    plot(x, theta, 'o-', 'LineWidth', 1.5);
%    grid on;
%    xlabel('x');
%    ylabel('\theta(x)');
%    title('Nodal rotation');

end

% ===== element routine =====
function Ke = beamHermiteKe(EI, le)
%BEAMHERMITEKE  4x4 Euler–Bernoulli beam element stiffness (Hermite cubic).
% DOFs: [w1, th1, w2, th2]
    c = EI / le^3;
    Ke = c * [ ...
        12,      6*le,   -12,      6*le; ...
        6*le,  4*le^2,  -6*le,   2*le^2; ...
       -12,     -6*le,    12,     -6*le; ...
        6*le,  2*le^2,  -6*le,   4*le^2  ...
    ];
end






