function solver(f::Function,
                hdiag::Vector{T},
                nroots::Int64,
                matdim::Int64;
                tol=1e-4,
                blocksize=nroots+5,
                maxvec=4*blocksize,
                niter=100,
                verbose=true
                ) where T <: AllowedFloat

    # Check on the input
    checkinp(nroots, blocksize, maxvec)

    # Davidson cache
    cache = DavidsonCache{T}(f, hdiag, nroots, matdim, blocksize,
                                   maxvec, tol, niter)

    # Construct the guess vectors
    guessvec(cache)

    # Run the generalised Davidson algorithm
    run_gendav(cache, verbose::Bool)

    # Eigenvalues
    values = cache.rho[1:nroots]
    
    # Get the eigenvectors
    vectors = Matrix{T}(undef, matdim, nroots)
    eigenvectors!(vectors, cache)
    
    return vectors, values
    
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

function run_gendav(cache::Cache, verbose::Bool)

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
        
        # Compute the new elements in the subspace matrix
        subspace_matrix(cache)

        # Compute the eigenpairs of the subspace matrix
        subspace_diag(cache)

        # Compute the residual vectors. Note that these will be stored
        # in the bvec array and subsequently transformed in place
        # to obtain the correction vectors, and then the new
        # subspace vectors
        residual_vectors(cache)

        # Print the report for this iteration
        if verbose print_report(k, cache) end

        # Stop here if all the roots are converged
        if all(i -> i == 1, cache.iconv[1:cache.nroots]) break end

        # Compute the correction vectors
        correction_vectors(cache)

        # Compute the new subspace vectors
        subspace_vectors(cache)

        # Subspace collapse?
        if cache.currdim + cache.nnew * 2 > cache.maxvec
            subspace_collapse(cache)
        else
            cache.currdim = cache.currdim + cache.nnew
        end

        # If we are here and this is the last iteration, then
        # we failed to converge all roots
        if k == cache.niter @warn "Not all roots converged" end
        
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

function subspace_matrix(cache::DavidsonCache{T}
                         ) where T <: AbstractFloat

    # Dimensions
    @unpack currdim, nnew, zero, one = cache

    # Work array
    bsigma = reshape(view(cache.work3, 1:currdim*nnew),
                     (currdim, nnew))
    fill!(bsigma, 0.0)

    # b^T sigma matrix product
    i1 = currdim - nnew + 1
    i2 = currdim
    b = view(cache.bvec, :, 1:i2)
    σ = view(cache.sigvec, :, i1:i2)
    BLAS.gemm!('T', 'N', one, b, σ, zero, bsigma)

    # Fill in the Gmat array
    for i in 1:cache.currdim
        for j in 1:cache.nnew
            j1 = cache.currdim - cache.nnew + j
            cache.Gmat[i,j1] = bsigma[i,j]
            cache.Gmat[j1,i] = cache.Gmat[i,j1]
        end
    end

end

function subspace_diag(cache::DavidsonCache{T}) where T <: AllowedFloat

    currdim = cache.currdim
    
    JOBZ = "V"
    UPLO = "L"
    N = currdim
    A = reshape(view(cache.alpha, 1:currdim^2), (currdim, currdim))
    LDA = currdim
    W = view(cache.rho, 1:currdim)
    LWORK = cache.lwork
    WORK = view(cache.evwork, 1:LWORK)
    INFO = cache.info

    # Fill in the subspace matrix
    G = cache.Gmat
    for j in 1:currdim
        for i in 1:currdim
            A[i,j] = G[i,j]
        end
    end

    # Call to ?syev
    syev!(JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO)
    
end



function residual_vectors(cache::Cache)

    @unpack matdim, currdim, blocksize, nnew, tol, minus_one, zero, one = cache

    # Subspace eigenvectors    
    alpha = reshape(view(cache.alpha, 1:currdim^2), (currdim, currdim))

    # Subspace eigenvalues
    rho = view(cache.rho, 1:currdim)
    
    # Working arrays
    alpha_bar = reshape(view(cache.work3, 1:currdim*blocksize),
                        (currdim, blocksize))
    
    #
    # Compute the residual vectors
    # r_k = Sum_i alpha_ik * (sigma_i - rho_k b_i),
    #
    # (alpha_bar)_ik = alpha_ik * rho_K
    for k in 1:blocksize
        for i in 1:currdim
            alpha_bar[i,k] = alpha[i,k] * rho[k]
        end
    end
    
    # sigma alpha
    σ = view(cache.sigvec, 1:matdim, 1:currdim)
    α = reshape(view(cache.alpha, 1:currdim*blocksize),
                (currdim, blocksize))
    work2 = cache.work2
    BLAS.gemm!('N', 'N', one, σ, α, zero, work2)
    
    # sigma alpha - b alpha_bar
    bvec = view(cache.bvec, 1:matdim, 1:currdim)
    BLAS.gemm!('N', 'N', minus_one, bvec, alpha_bar, one, work2)
        
    #
    # Save the residual vectors for the unconverged roots
    #
    # Initialisation
    ki=currdim + 1
    kf=currdim + nnew

    bvec = cache.bvec
    rnorm = cache.rnorm
    iconv = cache.iconv
    rho1 = cache.rho1
    
    # Loop over roots
    nnew = 0
    iconv .= 0

    for k in 1:blocksize

        # Residual norm
        w2 = view(work2, :, k:k)
        rnorm[k] = sqrt(dot(w2, w2))
        
        # Update the convergence information
        if rnorm[k] < tol iconv[k] = 1 end
    
        # Save the residual vector and corresponding eigenvalue if it
        # corresponds to an unconverged root
        if iconv[k] == 0

            nnew = nnew + 1

            k1 = ki-1+nnew
            b = view(bvec, :, k1:k1)
            w2 = view(work2, :, k:k)
            b = copy!(b, w2) 
            
            rho1[nnew] = rho[k]

        end

    end
    
end

function print_report(k::Int64, cache::Cache)

    # Table header
    if k == 1
        println("\n", "*"^43)
        println(" Iteration  Nvec   Max rnorm       Nconv")
        println("*"^43)
    end

    # Information for the current iteration

    #if (verbose) write(6,'(x,i4,7x,i4,3x,E13.7,3x,i4)') &
    #     k,currdim,maxval(rnorm(1:nstates)),sum(iconv(1:nstates))

    @unpack currdim, rnorm, iconv, nroots = cache
    maxres = maximum(rnorm[1:nroots])
    nconv = sum(iconv)

    println("$k $currdim $maxres $nconv")
    
end

function correction_vectors(cache::Cache)

    @unpack matdim, currdim, nnew, hdiag, rho1 = cache

    bvec = cache.bvec
    
    #
    # Diagonal preconditioned residue correction vectors
    #

    # Indices of the positions in the bvec array in which the
    # correction vectors will be stored
    ki = currdim + 1
    kf = currdim + nnew
    
    # Loop over correction vectors
    k1 = 0
    for k in ki:kf

        k1 += 1
        
        # Loop over elements of the correction vector
        for i in 1:matdim
            bvec[i,k] = -bvec[i,k] / (hdiag[i]-rho1[k1])
        end
            
    end

end

function subspace_vectors(cache::Cache)

    @unpack T, matdim, currdim, nnew, zero, one = cache

    bvec = cache.bvec
    sigvec = cache.sigvec
    work2 = cache.work2
    
    # Indices of the positions in the bvec array in which the new
    # subspace vectors will be stored
    ki = currdim + 1
    kf = currdim + nnew
    
    #
    # New orthonormalisation of the correction vectors
    #
    # Performed in two steps:
    #
    # (1) Gram-Schmidt orthogonalisation against the previous subspace
    #     vectors
    #
    # (2) Symmetric orthogonalisation within the space spanned by the
    #     intermediately orthogonalised correction vectors from (1)
    #

    # Overlaps between the previous subspace vectors and the correction
    # vectors
    Smat = reshape(view(cache.work3, 1:currdim*nnew), (currdim, nnew))
    
    bprev = view(bvec, :, 1:currdim)
    bnew = view(bvec, :, ki:kf)
    BLAS.gemm!('T', 'N', one, bprev, bnew, zero, Smat)

    # GS orthogonalisation of the correction vectors against the previous
    # subspace vectors
    k1 = 0
    for k in ki:kf
        k1 += 1
        for i in 1:currdim
            for j in 1:matdim
                bvec[j,k] -= Smat[i,k1] * bvec[j,i]
            end
            
        end
        bk = reshape(view(bvec, :, k:k), (matdim))
        len_bk = sqrt(dot(bk, bk))
        for j in 1:matdim
            bk[j] /= len_bk
        end
    end
    
    # Overlaps between the intermediately orthogonalised correction
    # vectors
    Smat = reshape(view(cache.work3, 1:nnew*nnew), (nnew, nnew))
    BLAS.gemm!('T', 'N', one, bnew, bnew, zero, Smat)

    # Inverse square root of the overlap matrix
    # Note that invsqrt_matrix will overwrite the Smat matrix
    Sinvsq = reshape(view(cache.work4, 1:nnew*nnew), (nnew, nnew))
    invsqrt_matrix!(Sinvsq, Smat, cache)
    
    # Symmetric orthogonalisation of the intermediately orthogonalised
    # correction vectors amongst themselves
    w2 = view(work2, :, 1:nnew)
    BLAS.gemm!('N', 'N', one, bnew, Sinvsq, zero, w2)
    copy!(bnew, w2)
        
end

function invsqrt_matrix!(Ainvsq::AbstractMatrix{T},
                         A::AbstractMatrix{T}, cache::Cache
                         ) where T <: AllowedFloat

    #
    # N.B. this will overwrite the contents of the input matrix A
    #

    # Dimension of the input matrix
    dim = size(A)[1]

    # Eigenpairs of the input matrix
    JOBZ = "V"
    UPLO = "L"
    N = dim
    LDA = dim
    W = view(cache.work1, 1:dim)
    LWORK = cache.lwork
    WORK = view(cache.evwork, 1:LWORK)
    INFO = cache.info
    syev!(JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO)

    # Inverse square root of the input matrix
    fill!(Ainvsq, 0.0)
    for k in 1:dim
        λinvsq = 1.0 / sqrt(abs(W[k]))
        for j in 1:dim
            for i in 1:dim
                Ainvsq[j,i] += A[i,k] * A[j,k] * λinvsq
            end
        end
    end
    
end
    
function subspace_collapse(cache::Cache)

    #
    # Compute the Ritz vectors. We will use the sigvec array as a
    # working array here to store the Ritz vectors
    #
    
    # Initialisation
    blocksize = cache.blocksize
    currdim = cache.currdim
    sigvec = cache.sigvec
    fill!(sigvec, 0.0)
    alpha = reshape(view(cache.alpha, 1:currdim^2), (currdim, currdim))
    bvec = cache.bvec
    
    # Compute the Ritz vectors
    @unpack one, zero = cache
    b = view(bvec, :, 1:currdim)
    α = view(alpha, 1:currdim, 1:blocksize)
    σ = view(sigvec, :, 1:blocksize)
    BLAS.gemm!('N', 'N', one, b, α, zero, σ)
            
    #
    # Collapse the subspace to be spanned by the lowest-lying Ritz vectors
    #

    # New subspace dimension
    cache.currdim = cache.blocksize
    cache.nnew = cache.blocksize

    # Save the Ritz vectors as the new subspace vectors
    b = view(bvec, :, 1:blocksize)
    σ = view(sigvec, :, 1:blocksize)
    copy!(b, σ)
    
    
end

function eigenvectors!(vectors::Matrix{T}, cache::Cache
                       ) where T <: AllowedFloat
    
    # Compute the Ritz vectors for the nroots lowest roots
    @unpack currdim, nroots, zero, one = cache
    alpha = reshape(view(cache.alpha, 1:currdim*nroots),
                    (currdim, nroots))
    bvec = view(cache.bvec, :, 1:currdim)

    BLAS.gemm!('N', 'N', one, bvec, alpha, zero, vectors)
    
end

