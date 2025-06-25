function solver(f::Function, hdiag::Vector{Float64}, nroots::Int64,
                matdim::Int64; tol=1e-4, blocksize=nroots+5,
                maxvec=4*blocksize, niter=100)
    
    # Check on the input
    checkinp(nroots, blocksize, maxvec)

    # Davidson cache
    cache = DavidsonCache(f, hdiag, nroots, matdim, blocksize, maxvec,
                          tol, niter)

    # Construct the guess vectors
    guessvec(cache)

    # Run the generalised Davidson algorithm
    run_gendav(cache)
    
end

function checkinp(nroots::Int64, blocksize::Int64, maxvec::Int64)

    # The block size must be greater than or equal to the number
    # of roots
    if blocksize < nroots
        @error "Block size is less than the no. roots"
    end

    # The maximum subspace dimension must be a multiple of the
    # block size
    if mod(maxvec, blocksize) != 0
        @error "The maximum subspace dimension must be a multiple" *
            " of the block size"
    end
    
end

function guessvec(cache::Cache)

    # For now, we will take the unit vectors (0,...,0,1,0,..,0)
    # as our guess vectors with the non-zero elements corresponding
    # to the lowest value diagonal matrix elements

    @unpack f, hdiag, nroots, matdim, blocksize, maxvec, tol, niter,
    bvec, sigvec, Gmat, alpha, rho, rho1, rnorm, work, work2, currdim,
    nconv, nnew, nsigma, iconv = cache
    
    # Sort the diagonal matrix elements
    ix = sortperm(hdiag)

    # Construct the guess vectors
    fill!(bvec, 0.0)
    for i in blocksize
        bvec[ix[i],i] = 1.0
    end
    
end

function run_gendav(cache::Cache)
    
    @unpack f, hdiag, nroots, matdim, blocksize, maxvec, tol, niter,
    bvec, sigvec, Gmat, alpha, rho, rho1, rnorm, work, work2, currdim,
    nconv, nnew, nsigma, iconv = cache

    #
    # Initialisation
    #
    # currdim: the current dimension of the subspace
    # nnew:    the no. new subspace vectors added in a given iteration
    # nconv:   the no. of converged roots
    # nsigma:  the total no. sigma vectors calculated
    currdim=blocksize
    nnew=blocksize
    nconv=0
    nsigma=0

    # Loop over iterations
    for k in 1:niter

        # Compute the σ-vectors
        sigma_vectors(cache)
        
        # Compute the new elements in the subspace Hamiltonian
        # matrix
        subspace_hamiltonian(cache)
        
    end
    
end

function sigma_vectors(cache::Cache)

    @unpack f, hdiag, nroots, matdim, blocksize, maxvec, tol, niter,
    bvec, sigvec, Gmat, alpha, rho, rho1, rnorm, work, work2, currdim,
    nconv, nnew, nsigma, iconv = cache
    
    # Update the total no. σ-vector calculations
    nsigma = nsigma + nnew
    
    # Indices of the first and last subspace vectors for which
    # σ-vectors are required
    ki = currdim - nnew + 1
    kf = currdim
    
    # Compute the σ-vectors
    f(bvec[:,ki:kf], sigvec[:,ki:kf])
    
end

function subspace_hamiltonian(cache::Cache)

    @unpack f, hdiag, nroots, matdim, blocksize, maxvec, tol, niter,
    bvec, sigvec, Gmat, alpha, rho, rho1, rnorm, work, work2, currdim,
    nconv, nnew, nsigma, iconv = cache
    
    # Work array
    bsigma = Matrix{Float64}(undef, currdim, nnew)

    # b^T sigma matrix product
    i1 = currdim - nnew + 1
    i2 = currdim

    
    
    exit()
    
end
