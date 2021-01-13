"""
    Computes mass and stress for the 10-bar truss structure

# Arguments
- `A::Vector{Float64}`: cross-sectional areas of each bar.
    See image in book for number order if needed.

# Returns
- `mass::Float64`: mass of the entire structure
- `stress::Vector{Float64}`: corresponding stress in each bar
"""
function truss(A)    
    
    P = 1e5  # applied loads
    Ls = 360  # length of sides
    Ld = sqrt(360^2 * 2)  # length of diagonals
    
    start = [5, 3, 6, 4, 4, 2, 5, 6, 3, 4]
    finish = [3, 1, 4, 2, 3, 1, 4, 3, 2, 1]
    phi = [0, 0, 0, 0, 90, 90, -45, 45, -45, 45]*pi/180
    L = [Ls, Ls, Ls, Ls, Ls, Ls, Ld, Ld, Ld, Ld]
    
    nbar = length(A)  # number of bars
    E = 1e7*ones(nbar)  # modulus of elasticity
    rho = 0.1*ones(nbar)  # material density
    
    Fx = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    Fy = [0.0, -P, 0.0, -P, 0.0, 0.0]
    # rigid = [0, 0, 0, 0, 1, 1]
    rigid = [false, false, false, false, true, true]
    
    n = length(Fx)  # number of nodes
    DOF = 2  # number of degrees of freedom
    
    # compute mass
    mass = sum(rho.*A.*L)
    
    # assemble global matricies
    K = zeros(DOF*n, DOF*n)
    S = zeros(nbar, DOF*n)
    
    for i = 1:nbar  # loop through each bar
    
        # compute the submatrices for the element
        Ksub, Ssub = bar(E[i], A[i], L[i], phi[i])
    
        # insert submatrix into global matrix
        idx = node2idx([start[i], finish[i]], DOF)  # pass in the starting and ending node number for this element
        K[idx, idx] += Ksub
        S[i, idx] = Ssub
    end
    
    # setup applied loads
    F = zeros(n*DOF)
    
    for i = 1:n
        idx = node2idx(i, DOF)
        F[idx[1]] = Fx[i]
        F[idx[2]] = Fy[i]
    end
    
    # setup boundary condition
    remove = node2idx(findall(rigid), DOF)
    keep = setdiff(1:n*DOF, remove)

    K = K[keep, keep]
    F = F[keep]
    S = S[:, keep]
    
    # solve for deflections
    d = K\F
    
    # compute stress
    stress = S*d

    return mass, stress
end
    

"""
    Compute the stiffness and stress matrix for one element

# Arguments
- `E::Float64`: modulus of elasticity
- `A::Float64`: cross-sectional area
- `L::Float64`: length of element
- `phi::Float64`: orientation of element

# Returns
- `K::Matrix{Float64}`: 4x4 stiffness matrix
- `S::Matrix{Float64}`: 1x4 stress matrix
"""
function bar(E, A, L, phi)
    
    # rename
    c = cos(phi)
    s = sin(phi)
    
    # stiffness matrix
    k0 = [c^2 c*s
          c*s s^2]
    
    K = E*A/L*[k0 -k0
               -k0 k0]
    
    # stress matrix
    S = E/L*[-c -s c s]

    return K, S
end

"""
Computes the appropriate indices in the global matrix for
the corresponding node numbers.  You pass in the number of the node
(either as a scalar or an array of locations), and the degrees of
freedom per node and it returns the corresponding indices in
the global matrices
"""
function node2idx(node, DOF)
    
    idx = Int64[]
    
    for i = 1:length(node)
        start = DOF*(node[i]-1) + 1
        finish = DOF*node[i]
    
        idx = [idx; start:finish]
    end
    return idx
end