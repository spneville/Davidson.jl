mutable struct DavidsonCache <: Cache

    # Matrix-vector multiplication function
    f::Function

    # Diagonal of the matrix whose eigenpairs are sought
    hdiag::Vector{Float64}
    
    # Number of roots to compute
    nroots::Int64

    # Dimension of the matrix whose eigenpairs are sought
    matdim::Int64

    # Block size
    blocksize::Int64

    # Maximum subspace dimension
    maxvec::Int64

    # Residual norm convergence threshold
    tol::Float64

    # Maximum number of iterations
    niter::Int64

    # Subspace vectors
    bvec::Matrix{Float64}
        
    # Sigma vectors
    sigvec::Matrix{Float64}
    
    # Subspace Hamiltonian matrix and eigenpairs
    Gmat::Matrix{Float64}
    alpha::Matrix{Float64}
    rho::Vector{Float64}
    rho1::Vector{Float64}
    
    # Residual norms
    rnorm::Vector{Float64}
    
    # Work arrays
    work::Vector{Float64}
    work2::Matrix{Float64}
    
    # Counters, etc.
    currdim::Int64
    nconv::Int64
    nnew::Int64
    nsigma::Int64
    iconv::Vector{Int64}
    
    # Inner constructor
    function DavidsonCache(f::Function, hdiag::Vector{Float64},
                           nroots::Int64, matdim::Int64,
                           blocksize::Int64, maxvec::Int64, tol::Float64,
                           niter::Int64)

        # Subspace vectors
        bvec = Matrix{Float64}(undef, matdim, maxvec)
        
        # Sigma vectors
        sigvec = Matrix{Float64}(undef, matdim, maxvec)
        
        # Subspace Hamiltonian matrix and eigenpairs
        Gmat = Matrix{Float64}(undef, maxvec,maxvec)
        alpha = Matrix{Float64}(undef, maxvec,maxvec)
        rho = Vector{Float64}(undef, maxvec)
        rho1 = Vector{Float64}(undef, maxvec)
        
        # Residual norms
        rnorm = Vector{Float64}(undef, blocksize)
        
        # Work arrays
        work = Vector{Float64}(undef, matdim)
        work2 = Matrix{Float64}(undef, matdim,blocksize)
        
        # Counters, etc.
        currdim = 0
        nconv = 0
        nnew = 0
        nsigma = 0
        iconv = Vector{Int64}(undef, blocksize)
        
        new(f, hdiag, nroots, matdim, blocksize, maxvec, tol, niter,
            bvec, sigvec, Gmat, alpha, rho, rho1, rnorm, work, work2,
            currdim, nconv, nnew, nsigma, iconv)
        
    end
        
end
