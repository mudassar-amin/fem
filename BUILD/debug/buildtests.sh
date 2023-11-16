#!/bin/bash

if [[ $# != 1 || $1 == *help ]]
then
  echo "usage: $0 regexp"
  echo "  Builds tests matching the regexp."
  echo "  The EIGEN_MAKE_ARGS environment variable allows to pass args to 'make'."
  echo "    For example, to launch 5 concurrent builds, use EIGEN_MAKE_ARGS='-j5'"
  exit 0
fi

TESTSLIST="rand
meta
maxsizevector
numext
sizeof
dynalloc
nomalloc
first_aligned
type_alias
nullary
mixingtypes
io
packetmath
vectorization_logic
basicstuff
constexpr
constructor
linearstructure
integer_types
unalignedcount
exceptions
redux
visitor
block
corners
symbolic_index
indexed_view
reshape
swap
resize
conservative_resize
product_small
product_large
product_extra
diagonalmatrices
skew_symmetric_matrix3
adjoint
diagonal
miscmatrices
commainitializer
smallvectors
mapped_matrix
mapstride
unaryviewstride
mapstaticmethods
array_cwise
array_for_matrix
array_replicate
array_reverse
ref
is_same_dense
triangular
selfadjoint
product_selfadjoint
product_symm
product_syrk
product_trmv
product_trmm
product_trsolve
product_mmtr
product_notemporary
stable_norm
permutationmatrices
bandmatrix
cholesky
lu
determinant
inverse
qr
qr_colpivoting
qr_fullpivoting
upperbidiagonalization
hessenberg
schur_real
schur_complex
eigensolver_selfadjoint
eigensolver_generic
eigensolver_complex
real_qz
eigensolver_generalized_real
jacobi
jacobisvd
bdcsvd
householder
geo_orthomethods
geo_quaternion
geo_eulerangles
geo_parametrizedline
geo_alignedbox
geo_hyperplane
geo_transformations
geo_homogeneous
stdvector
stdvector_overload
stdlist
stdlist_overload
stddeque
stddeque_overload
sparse_basic
sparse_block
sparse_vector
sparse_product
sparse_ref
sparse_solvers
sparse_permutations
simplicial_cholesky
conjugate_gradient
incomplete_cholesky
bicgstab
lscg
sparselu
sparseqr
umeyama
nesting_ops
nestbyvalue
zerosized
dontalign
evaluators
sizeoverflow
prec_inverse_4x4
vectorwiseop
special_numbers
rvalue_types
dense_storage
ctorleak
inplace_decomposition
half_float
bfloat16_float
array_of_string
num_dimensions
stl_iterators
blasutil
random_matrix
initializer_list_construction
diagonal_matrix_variadic_ctor
serializer
tuple_test
threads_eventcount
threads_runqueue
threads_non_blocking_thread_pool
fastmath
boostmultiprec
NonLinearOptimization
NumericalDiff
autodiff_scalar
autodiff
BVH
matrix_exponential
matrix_function
matrix_power
matrix_square_root
alignedvector3
FFT
EulerAngles
NNLS
sparse_extra
polynomialsolver
polynomialutils
splines
gmres
dgmres
minres
idrs
bicgstabl
idrstabl
levenberg_marquardt
kronecker_product
bessel_functions
special_functions
special_packetmath
cxx11_tensor_argmax
cxx11_tensor_assign
cxx11_tensor_block_access
cxx11_tensor_block_eval
cxx11_tensor_block_io
cxx11_tensor_broadcasting
cxx11_tensor_casts
cxx11_tensor_chipping
cxx11_tensor_comparisons
cxx11_tensor_concatenation
cxx11_tensor_const
cxx11_tensor_contraction
cxx11_tensor_convolution
cxx11_tensor_custom_index
cxx11_tensor_custom_op
cxx11_tensor_dimension
cxx11_tensor_empty
cxx11_tensor_executor
cxx11_tensor_expr
cxx11_tensor_fft
cxx11_tensor_fixed_size
cxx11_tensor_forced_eval
cxx11_tensor_generator
cxx11_tensor_ifft
cxx11_tensor_image_patch
cxx11_tensor_index_list
cxx11_tensor_inflation
cxx11_tensor_intdiv
cxx11_tensor_io
cxx11_tensor_layout_swap
cxx11_tensor_lvalue
cxx11_tensor_map
cxx11_tensor_math
cxx11_tensor_mixed_indices
cxx11_tensor_morphing
cxx11_tensor_move
cxx11_tensor_notification
cxx11_tensor_of_complex
cxx11_tensor_of_const_values
cxx11_tensor_of_strings
cxx11_tensor_padding
cxx11_tensor_patch
cxx11_tensor_random
cxx11_tensor_reduction
cxx11_tensor_ref
cxx11_tensor_roundings
cxx11_tensor_scan
cxx11_tensor_shuffling
cxx11_tensor_simple
cxx11_tensor_striding
cxx11_tensor_sugar
cxx11_tensor_thread_local
cxx11_tensor_thread_pool
cxx11_tensor_trace
cxx11_tensor_volume_patch
cxx11_tensor_uint128
"
targets_to_make=$(echo "$TESTSLIST" | grep -E "$1" | xargs echo)

if [ -n "${EIGEN_MAKE_ARGS:+x}" ]
then
  C:/mingw64/bin/ninja.exe $targets_to_make ${EIGEN_MAKE_ARGS}
else
  C:/mingw64/bin/ninja.exe $targets_to_make 
fi
exit $?

