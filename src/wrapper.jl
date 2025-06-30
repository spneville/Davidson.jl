function syev!(JOBZ::String,
               UPLO::String,
               N::Int64,
               A::AbstractMatrix{Float64},
               LDA::Int64,
               W::AbstractVector{Float64},
               WORK::AbstractVector{Float64},
               LWORK::Int32,
               INFO::Int32)


    ccall((:dsyev_, "libblastrampoline"),
          Cvoid,
          (Cstring,      # JOBZ
           Cstring,      # UPLO
           Ref{Int32},   # N
           Ptr{Float64}, # A
           Ref{Int32},   # LDA
           Ptr{Float64}, # W
           Ptr{Float64}, # WORK
           Ref{Int32},   # LWORK
           Ref{Int32}),  #INFO
          JOBZ,
          UPLO,
          N,
          A,
          LDA,
          W,
          WORK,
          LWORK,
          INFO)
    
    
end

function syev!(JOBZ::String,
               UPLO::String,
               N::Int64,
               A::AbstractMatrix{Float32},
               LDA::Int64,
               W::AbstractVector{Float32},
               WORK::AbstractVector{Float32},
               LWORK::Int32,
               INFO::Int32)


    ccall((:ssyev_, "libblastrampoline"),
          Cvoid,
          (Cstring,      # JOBZ
           Cstring,      # UPLO
           Ref{Int32},   # N
           Ptr{Float64}, # A
           Ref{Int32},   # LDA
           Ptr{Float64}, # W
           Ptr{Float64}, # WORK
           Ref{Int32},   # LWORK
           Ref{Int32}),  #INFO
          JOBZ,
          UPLO,
          N,
          A,
          LDA,
          W,
          WORK,
          LWORK,
          INFO)
    
    
end

