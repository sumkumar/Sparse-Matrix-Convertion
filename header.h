typedef struct {
	long row, col, val;
}coo_tuple;

typedef struct {
	coo_tuple *index;
	long row_size, col_size, NZ_size;
}COO;

typedef struct {
	long *row_ptr, *col_ind, *val, row_size, col_size, NZ_size;
}CSR;

typedef struct {
	long *row_ind, *col_ptr, *val, row_size, col_size, NZ_size;
}CSC;

typedef struct COO_LL {
    coo_tuple coo_entry;
    struct COO_LL* next;
}COO_LL;

COO create_random_sparse_matrix(long m, long n, long elems, long neg_nums);

CSR create_csr_from_coo(COO* coo_mat);

CSC create_csc_from_coo(COO* coo_mat);

void csc_row_transform(CSC* csc_mat, long r1, long k1, long r2, long k2);

void print_coo_matrix(COO (*matrix));

void print_csc_matrix(CSC (*matrix));

void print_csr_matrix(CSR (*matrix));

COO csc_mul_csr(CSC* csc_mat, CSR* csr_mat);

COO csr_mul_csc(CSR* csr_mat, CSC* csc_mat);

COO create_coo_matrix(long row_size, long column_size, long NZ_size);

void coo_add_tuple(COO* matrix, long row, long col, long data, long index);
