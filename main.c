#include<stdio.h>
#include "header.h"

COO create_random_coo_matrix(){
    long i, m, n, row, column, data;
	double den_NZ, den_NN, elems, neg_nums_per, neg_nums;
	COO coo_mat;
	printf("Enter row size : \n");
	scanf("%ld", &m);
	if(m < 0){
        printf("Invalid Entry\n");
		coo_mat.NZ_size = -1;
        return coo_mat;
    }
	printf("Enter column size : \n");
	scanf("%ld", &n);
	if(n < 0){
        printf("Invalid Entry\n");
		coo_mat.NZ_size = -1;
        return coo_mat;
    }
	printf("Enter density of non-zero entries : \n");
	scanf("%lf", &den_NZ);
	printf("Enter percentage of negative numbers : \n");
	scanf("%lf", &neg_nums_per);
	den_NZ/=100;
	neg_nums_per/=100;
	elems = den_NZ*m*n;
	neg_nums = neg_nums_per*m*n;
	if(neg_nums > elems){
        printf("Negative numbers cannot be more than required numbers\n");
		coo_mat.NZ_size = -1;
        return coo_mat;
	}
	if(!((den_NZ >= 0 && den_NZ <= 1) && (neg_nums_per >=0 && neg_nums_per <= 1))){
        printf("Please give proper percentage values\n");
		coo_mat.NZ_size = -1;
        return coo_mat;
	}
    coo_mat = create_random_sparse_matrix(m, n, (int)elems, (int)neg_nums);
    /*for(i=0;i<coo_mat.NZ_size;i++){
        printf("row: %ld col: %ld val: %ld\n", coo_mat.index[i].row, coo_mat.index[i].col, coo_mat.index[i].val);
    }
    print_coo_matrix(&coo_test);*/
    return coo_mat;
}

COO get_coo_input(){
    long NZ_size, row_size, col_size, i, row, col, val;
    printf("Enter no. of Non-Zero integers: \n");
    scanf("%ld", &NZ_size);
    printf("Enter row size: ");
    scanf("%ld", &row_size);
    printf("Enter column size: ");
    scanf("%ld", &col_size);
    COO coo_mat = create_coo_matrix(row_size, col_size, NZ_size);
    for(i=0;i<NZ_size;i++){
        printf("Enter row no. of integer %ld : \n",i+1);
        scanf("%ld", &row);
        if(row > row_size || row < 0){
            printf("Invalid Entry\n");
			coo_mat.NZ_size = -1;
            return coo_mat;
        }
        printf("Enter column no. of integer %ld : \n",i+1);
        scanf("%ld", &col);
        if(col > col_size || col < 0){
            printf("Invalid Entry\n");
			coo_mat.NZ_size = -1;
            return coo_mat;
        }
        printf("Enter value of integer %ld : \n",i+1);
        scanf("%ld", &val);
        coo_add_tuple(&coo_mat, row, col, val, i);
    }
    return coo_mat;
}

CSR get_csr_from_coo(){
    COO coo_test = get_coo_input();
    return create_csr_from_coo(&coo_test);
}

CSC get_csc_from_coo(){
    COO coo_test = get_coo_input();
    return create_csc_from_coo(&coo_test);
}

void csc_row_transformation(CSC* csc_mat){
    long r1, r2, k1, k2;
    printf("Enter value of k1 : \n");
    scanf("%ld", &k1);
    printf("Enter value of r1 : \n");
    scanf("%ld", &r1);
    if(r1 > csc_mat->row_size || r1 < 0){
        printf("Invalid Entry\n");
        return ;
    }
    printf("Enter value of k2 : \n");
    scanf("%ld", &k2);
    printf("Enter value of r2 : \n");
    scanf("%ld", &r2);
    if(r2 > csc_mat->row_size || r2 < 0){
        printf("Invalid Entry\n");
        return ;
    }
    csc_row_transform(csc_mat, r1, k1, r2, k2);
}

void main(){
    printf("Enter information for Random COO matrix\n");
    COO coo_mat = create_random_coo_matrix();
	if(coo_mat.NZ_size == -1)
		return ;
    printf("Random COO matrix\n");
    print_coo_matrix(&coo_mat);
    CSC csc_mat = create_csc_from_coo(&coo_mat);
    printf("CSC matrix created from above Random COO matrix\n");
    print_csc_matrix(&csc_mat);
    printf("Enter CSC row transformation information\n");
    csc_row_transformation(&csc_mat);
    printf("CSC matrix after row transformation\n");
    print_csc_matrix(&csc_mat);
    printf("Enter information for Random COO matrix\n");
    COO coo_mat2 = create_random_coo_matrix();
	if(coo_mat.NZ_size == -1)
		return ;
    printf("Another Random COO matrix\n");
    print_coo_matrix(&coo_mat2);
    CSR csr_mat = create_csr_from_coo(&coo_mat2);
    printf("CSR matrix created from above Random COO matrix\n");
    print_csr_matrix(&csr_mat);
    printf("Multiplying CSC matrix with CSR matrix\n");
    COO coo_mul_result = csc_mul_csr(&csc_mat, &csr_mat);
    if(&coo_mul_result != NULL){
        printf("Resultant matrix after Multiplying CSC matrix with CSR matrix\n");
        print_coo_matrix(&coo_mul_result);
    }
    return ;
}
