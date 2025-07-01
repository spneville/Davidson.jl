mutable struct DavidsonCache{T} <: Cache
    
    # Type parameter
    T::Type
    
    # Matrix-vector multiplication function
    f::Function

    # Diagonal of the matrix whose eigenpairs are sought
    hdiag::Vector{T}
    
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
    bvec::Matrix{T}
        
    # Sigma vectors
    sigvec::Matrix{T}
    
    # Subspace matrix and eigenpairs
    Gmat::Matrix{T}
    alpha::Vector{T}
    rho::Vector{T}
    rho1::Vector{T}
    
    # Residual norms
    rnorm::Vector{Float64}
    
    # Work arrays
    work1::Vector{T}
    work2::Matrix{T}
    work3::Vector{T}
    work4::Vector{T}

    # Counters, etc.
    currdim::Int64
    nconv::Int64
    nnew::Int64
    nsigma::Int64
    iconv::Vector{Int64}

    # LAPACK ?syev and ?heev work arrays
    lwork::Int32
    evwork::Vector{T}
    info::Int32

    # -1.0, 0.0, and +1.0
    minus_one::T
    zero::T
    one::T
    
    # Inner constructor
    function DavidsonCache{T}(f::Function,
                              hdiag::Vector{T},
                              nroots::Int64,
                              matdim::Int64,
                              blocksize::Int64,
                              maxvec::Int64,
                              tol::Float64,
                              niter::Int64
                              ) where T <: AllowedFloat

        # Subspace vectors
        bvec = Matrix{T}(undef, matdim, maxvec)
        
        # Sigma vectors
        sigvec = Matrix{T}(undef, matdim, maxvec)
        
        # Subspace matrix and eigenpairs
        Gmat = Matrix{T}(undef, maxvec, maxvec)
        alpha = Vector{T}(undef, maxvec*maxvec)
        rho = Vector{T}(undef, maxvec)
        rho1 = Vector{T}(undef, maxvec)
        
        # Residual norms
        rnorm = Vector{Float64}(undef, blocksize)
        
        # Work arrays
        work1 = Vector{T}(undef, matdim)
        work2 = Matrix{T}(undef, matdim,blocksize)
        work3 = Vector{T}(undef, maxvec*blocksize)
        work4 = Vector{T}(undef, blocksize*blocksize)
        
        # Counters, etc.
        currdim = 0
        nconv = 0
        nnew = 0
        nsigma = 0
        iconv = Vector{Int64}(undef, blocksize)

        # LAPACK ?syev and ?heev work arrays
        lwork::Int32 = 3 * maxvec
        evwork = Vector{T}(undef, lwork)
        info::Int32 = 0

        # -1.0, 0.0, and +1.0
        minus_one::T = -1.0
        zero::T = 0.0
        one::T = 1.0
        
        new{T}(T, f, hdiag, nroots, matdim, blocksize, maxvec, tol, niter,
               bvec, sigvec, Gmat, alpha, rho, rho1, rnorm, work1, work2,
               work3, work4, currdim, nconv, nnew, nsigma, iconv,
               lwork, evwork, info, minus_one, zero, one)
        
    end
        
end
