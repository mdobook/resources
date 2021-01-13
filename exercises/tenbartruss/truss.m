function [mass, stress] = truss(A)
%{
Computes mass and stress for the 10-bar truss structure

Parameters
----------
A : array of length 10 (row or column)
    cross-sectional areas of each bar.
    See image in book for number order if needed.

Outputs
-------
mass : float
    mass of the entire structure
stress : array of length 10
    corresponding stress of each bar
%}


P = 1e5;  % applied loads
Ls = 360;  % length of sides
Ld = sqrt(360^2 * 2);  % length of diagonals

start = [5, 3, 6, 4, 4, 2, 5, 6, 3, 4];
finish = [3, 1, 4, 2, 3, 1, 4, 3, 2, 1];
phi = [0, 0, 0, 0, 90, 90, -45, 45, -45, 45]*pi/180;
L = [Ls, Ls, Ls, Ls, Ls, Ls, Ld, Ld, Ld, Ld];

nbar = length(A);
E = 1e7*ones(1, nbar);  % modulus of elasticity
rho = 0.1*ones(1, nbar);  % material density

Fx = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
Fy = [0.0, -P, 0.0, -P, 0.0, 0.0];
rigidx = [0, 0, 0, 0, 1, 1];
rigidy = [0, 0, 0, 0, 1, 1];


n = length(Fx);  % number of nodes
DOF = 2;  % number of degrees of freedom
nbar = length(A);  % number of bars

% compute mass
if iscolumn(A)
    A = A';
end
mass = sum(rho.*A.*L);

% assemble global matricies
K = zeros(DOF*n, DOF*n);
S = zeros(nbar, DOF*n);

for i = 1:nbar  % loop through each bar

    % compute the submatrices for the element
    [Ksub, Ssub] = bar(E(i), A(i), L(i), phi(i));

    % insert submatrix into global matrix
    idx = node2idx([start(i), finish(i)], DOF);  % pass in the starting and ending node number for this element
    K(idx, idx) = K(idx, idx) + Ksub;
    S(i, idx) = Ssub;
end

% setup applied loads
F = zeros(n*DOF, 1);

for i = 1:n
    idx = node2idx(i, DOF);
    F(idx(1)) = Fx(i);
    F(idx(2)) = Fy(i);
end

% setup boundary condition
idxx = find(rigidx);
idxy = find(rigidy);
removex = node2idx(idxx, DOF);
removey = node2idx(idxy, DOF);
removex = removex(1:2:end);  % take only x components
removey = removey(2:2:end);  % take only y components
remove = [removex, removey];  % concatenate

K(remove, :) = [];
K(:, remove) = [];
F(remove) = [];
S(:, remove) = [];


% solve for deflections
d = K\F;

% compute stress
stress = S*d;

end


% -------------


function [K, S] = bar(E, A, L, phi)
%{
Compute the stiffness and stress matrix for one element

Parameters
----------
E : float
    modulus of elasticity
A : float
    cross-sectional area
L : float
    length of element
phi : float
    orientation of element

Outputs
-------
K : 4 x 4 ndarray
    stiffness matrix
S : 1 x 4 ndarray
    stress matrix
%}


% rename
c = cos(phi);
s = sin(phi);

% stiffness matrix
k0 = [c^2 c*s;
      c*s s^2];

K = E*A/L*[k0 -k0;
           -k0 k0];

% stress matrix
S = E/L*[-c -s c s];
end


% --------


function idx = node2idx(node, DOF)
%{
Computes the appropriate indices in the global matrix for
the corresponding node numbers.  You pass in the number of the node
(either as a scalar or an array of locations), and the degrees of
freedom per node and it returns the corresponding indices in
the global matrices
%}

idx = [];

for i = 1:length(node)

    start = DOF*(node(i)-1) + 1;
    finish = DOF*node(i);

    idx = [idx start:finish];

end

end
