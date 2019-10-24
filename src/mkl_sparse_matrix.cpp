#include "sparse_matrix.hpp"

string MKL_SparseType::matrixType(){ std::cout<<"CSR Matrix from MKL Library."; };
void MKL_SparseType::ConvertFromCOO(vector<int> &rows,vector<int> &cols,vector<complex<double> > &vals)
{
	vector<int> rows_; rows_.swap(rows); 
	vector<int> cols_; cols_.swap(cols); 
	vector<complex<double> > vals_; vals_.swap(vals); 

	sparse_matrix_t newMatrix;
	assert( mkl_sparse_z_create_coo (&newMatrix, SPARSE_INDEX_BASE_ZERO, numRows(),numCols(),rows_.size(), &rows_[0],&cols_[0],&vals_[0])== SPARSE_STATUS_SUCCESS );
	assert( mkl_sparse_convert_csr (newMatrix, SPARSE_OPERATION_NON_TRANSPOSE, &Matrix)== SPARSE_STATUS_SUCCESS );
	assert( mkl_sparse_destroy(newMatrix)== SPARSE_STATUS_SUCCESS );

	descr.type = SPARSE_MATRIX_TYPE_HERMITIAN;
	descr.mode = SPARSE_FILL_MODE_UPPER;
	descr.diag = SPARSE_DIAG_NON_UNIT;

};
	

void MKL_SparseType::Multiply(const complex<double>* a,const complex<double>* x,const complex<double>*  b, complex<double>* y)
{
	assert( mkl_sparse_z_mv(SPARSE_OPERATION_NON_TRANSPOSE, *a, Matrix ,descr, x, *b,y) == SPARSE_STATUS_SUCCESS );
};

void MKL_SparseType::Optimize()
{

//	std::cout<<*this->rows_[1]<<" "<<*this->rows_[2]<<std::endl;
/*		mkl_sparse_set_mv_hint(A,SPARSE_OPERATION_NON_TRANSPOSE,descr, (MKL_INT)numMul );
		mkl_sparse_set_dotmv_hint(A,SPARSE_OPERATION_NON_TRANSPOSE,descr, (MKL_INT)numMul );
		mkl_sparse_optimize (A);*/

};

