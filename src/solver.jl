function solver(f::Function, nroots::Int64, matdim::Int64;
                tol=1e-4, blocksize=nroots+5, maxvec=4*blocksize,
                niter=100)
    
    # Check on the input
    checkinp(nroots, blocksize, maxvec)

    # Davidson cache
    cache = DavidsonCache(f, nroots, matdim, blocksize, maxvec, tol,
                          niter)

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
