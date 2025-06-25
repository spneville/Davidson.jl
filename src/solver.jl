function solver(f::Function, hdiag::Vector{Float64}, nroots::Int64,
                matdim::Int64; tol=1e-4, blocksize=nroots+5,
                maxvec=4*blocksize, niter=100)
    
    # Check on the input
    checkinp(nroots, blocksize, maxvec)

    # Davidson cache
    cache = DavidsonCache{Float64}(f, hdiag, nroots, matdim, blocksize,
                                   maxvec, tol, niter)

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

    #
    # For now, we will take the unit vectors (0,...,0,1,0,..,0)
    # as our guess vectors with the non-zero elements corresponding
    # to the lowest value diagonal matrix elements
    #
    
    # Sort the diagonal matrix elements
    ix = sortperm(cache.hdiag)

    # Construct the guess vectors
    fill!(cache.bvec, 0.0)
    for i in 1:cache.blocksize
        cache.bvec[ix[i],i] = 1.0
    end
    
end

function run_gendav(cache::Cache)
    
    #
    # Initialisation
    #
    # currdim: the current dimension of the subspace
    # nnew:    the no. new subspace vectors added in a given iteration
    # nconv:   the no. of converged roots
    # nsigma:  the total no. sigma vectors calculated

    cache.currdim = cache.blocksize
    cache.nnew = cache.blocksize
    cache.nconv = 0
    cache.nsigma = 0

    # Loop over iterations
    for k in 1:cache.niter

        # Compute the σ-vectors
        sigma_vectors(cache)
        
        # Compute the new elements in the subspace Hamiltonian
        # matrix
        subspace_hamiltonian(cache)

        # Compute the eigenpairs of the subspace Hamiltonian
        subspace_diag(cache)
        
    end
    
end

function sigma_vectors(cache::Cache)

    # Update the total no. σ-vector calculations
    cache.nsigma = cache.nsigma + cache.nnew
    
    # Indices of the first and last subspace vectors for which
    # σ-vectors are required
    ki = cache.currdim - cache.nnew + 1
    kf = cache.currdim
    
    # Compute the σ-vectors
    b = view(cache.bvec, :, ki:kf)
    σ = view(cache.sigvec, :, ki:kf)
    cache.f(b, σ)

end

function subspace_hamiltonian(cache::Cache)

    # Work array
    bsigma = Matrix{Float64}(undef, cache.currdim, cache.nnew)
    fill!(bsigma, 0.0)
    
    # b^T sigma matrix product
    i1 = cache.currdim - cache.nnew + 1
    i2 = cache.currdim
    BLAS.gemm!('T', 'N', 1.0, cache.bvec[:,1:i2], cache.sigvec[:,i1:i2],
               0.0, bsigma)

    # Fill in the Gmat array
    for i in 1:cache.currdim
        for j in 1:cache.nnew
            j1 = cache.currdim - cache.nnew + j
            cache.Gmat[i,j1] = bsigma[i,j]
            cache.Gmat[j1,i] = cache.Gmat[i,j1]
        end
    end

end

function subspace_diag(cache::Cache)

    #
    # This needs to be replaced with a BLAS call using
    # pre-allocated work arrays. Something along the lines
    # of what is done in FastLapackInterface.jl, for example
    #
    F = eigen(cache.Gmat[1:cache.currdim,1:cache.currdim])

    exit()
    
end

