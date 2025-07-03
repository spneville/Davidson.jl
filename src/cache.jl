mutable struct DavidsonCache{T, R} <: Cache
    
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
    bvec::Vector{T}
    
    # Sigma vectors
    sigvec::Vector{T}

    # Subspace matrix and eigenpairs
    Gmat::Vector{T}
    alpha::Vector{T}
    rho::Vector{R}
    rho1::Vector{R}
    
    # Residual norms
    rnorm::Vector{Float64}
    
    # Work arrays
    work1::Vector{R}
    work2::Vector{T}
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
    revwork::Vector{R}
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
                              ) where T <: AllowedTypes

        # Real type
        R = T <: Allowed64 ? Float64 : Float32

        # Subspace vectors
        bvec = Vector{T}(undef, matdim*maxvec)
        
        # Sigma vectors
        sigvec = Vector{T}(undef, matdim*maxvec)

        # Subspace matrix and eigenpairs
        Gmat = Vector{T}(undef, maxvec*maxvec)
        alpha = Vector{T}(undef, maxvec*maxvec)
        rho = Vector{R}(undef, maxvec)
        rho1 = Vector{R}(undef, maxvec)
        
        # Residual norms
        rnorm = Vector{Float64}(undef, blocksize)
        
        # Work arrays
        work1 = Vector{R}(undef, matdim)
        work2 = Vector{T}(undef, matdim*blocksize)
        work3 = Vector{T}(undef, maxvec*blocksize)
        work4 = Vector{T}(undef, blocksize*blocksize)
        
        # Counters, etc.
        currdim = 0
        nconv = 0
        nnew = 0
        nsigma = 0
        iconv = Vector{Int64}(undef, blocksize)

        # LAPACK ?syev and ?heev work arrays
        # We will use the same dimension for both the work and rwork
        # arrays
        lwork::Int32 = 3 * maxvec
        evwork = Vector{T}(undef, lwork)
        revwork = Vector{R}(undef, lwork)
        info::Int32 = 0
        
        # -1.0, 0.0, and +1.0
        minus_one::T = -1.0
        zero::T = 0.0
        one::T = 1.0
        
        new{T, R}(T, f, hdiag, nroots, matdim, blocksize, maxvec, tol,
                  niter, bvec, sigvec, Gmat, alpha, rho, rho1, rnorm,
                  work1, work2, work3, work4, currdim, nconv, nnew,
                  nsigma, iconv, lwork, evwork, revwork, info, minus_one,
                  zero, one)
        
    end

end

function workarrays(T::DataType, matdim::Int64, blocksize::Int64,
                    maxvec::Int64)

    @assert T <: AllowedTypes
    
    if T <: Union{Float32, ComplexF32}
        R = Float32
    else
        R = Float64
    end
    
    Tdim = Tworksize(matdim, blocksize, maxvec)
    Twork = Vector{T}(undef, Tdim)
    
    Rdim = Rworksize(matdim, blocksize, maxvec)
    Rwork = Vector{R}(undef, Rdim)
    
    return Twork, Rwork
    
end

function Tworksize(matdim::Int64, blocksize::Int64, maxvec::Int64)

    dim = 0
    
    # Subspace vectors
    dim += matdim * maxvec
        
    # Sigma vectors
    dim += matdim * maxvec
    
    # Subspace matrix
    dim += maxvec * maxvec

    # Subspace eigenvectors
    dim += maxvec * maxvec

    # Residual norms
    dim += blocksize
    
    # Work2 array
    dim += matdim * blocksize

    # Work3 array
    dim += maxvec * blocksize

    # Work4 array
    dim += blocksize * blocksize
    
    # LAPACK ?syev and ?heev work arrays
    dim += 3 * maxvec

    return dim
    
end

function Rworksize(matdim::Int64, blocksize::Int64, maxvec::Int64)

    dim = 0
    
    # Subspace eigenvalues plus a working copy
    dim += maxvec
    dim += maxvec
    
    # Work1 array
    dim += matdim

    # LAPACK ?syev and ?heev work arrays
    dim += 3 * maxvec

    return dim
    
end
