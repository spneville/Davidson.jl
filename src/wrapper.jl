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

function heev!(JOBZ::String,
               UPLO::String,
               N::Int64,
               A::AbstractMatrix{ComplexF64},
               LDA::Int64,
               W::AbstractVector{Float64},
               WORK::AbstractVector{ComplexF64},
               LWORK::Int32,
               RWORK::AbstractVector{Float64},
               INFO::Int32)

    ccall((:zheev_, "libblastrampoline"),
          Cvoid,
          (Cstring,         # JOBZ
           Cstring,         # UPLO
           Ref{Int32},      # N
           Ptr{ComplexF64}, # A
           Ref{Int32},      # LDA
           Ptr{Float64},    # W
           Ptr{ComplexF64}, # WORK
           Ref{Int32},      # LWORK
           Ptr{Float64},    # RWORK
           Ref{Int32}),     #INFO
          JOBZ,
          UPLO,
          N,
          A,
          LDA,
          W,
          WORK,
          LWORK,
          RWORK,
          INFO)
    
end

function heev!(JOBZ::String,
               UPLO::String,
               N::Int64,
               A::AbstractMatrix{ComplexF32},
               LDA::Int64,
               W::AbstractVector{Float32},
               WORK::AbstractVector{ComplexF32},
               LWORK::Int32,
               RWORK::AbstractVector{Float32},
               INFO::Int32)

    ccall((:cheev_, "libblastrampoline"),
          Cvoid,
          (Cstring,         # JOBZ
           Cstring,         # UPLO
           Ref{Int32},      # N
           Ptr{ComplexF32}, # A
           Ref{Int32},      # LDA
           Ptr{Float32},    # W
           Ptr{ComplexF32}, # WORK
           Ref{Int32},      # LWORK
           Ptr{Float32},    # RWORK
           Ref{Int32}),     #INFO
          JOBZ,
          UPLO,
          N,
          A,
          LDA,
          W,
          WORK,
          LWORK,
          RWORK,
          INFO)
    
end
