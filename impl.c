#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include "header.h"

long print_mat[20][20];


CSR create_csr(long NZ_size, long row_size, long col_size){
	long* col_ind = (long*) malloc(NZ_size * sizeof(long));
	long* val = (long*) malloc(NZ_size * sizeof(long));
	long* row_ptr = (long*) malloc((row_size + 1) * sizeof(long));
	CSR new_csr = {row_ptr, col_ind, val, row_size, col_size, NZ_size};
	return new_csr;
}


CSC create_csc(long NZ_size, long row_size, long col_size){
	long* row_ind = (long*) malloc(NZ_size * sizeof(long));
	long* val = (long*) malloc(NZ_size * sizeof(long));
	long* col_ptr = (long*) malloc((col_size + 1) * sizeof(long));
	CSC new_csc = {row_ind, col_ptr, val, row_size, col_size, NZ_size};
	return new_csc;
}

COO create_coo_matrix(long row_size, long column_size, long NZ_size){
	COO mat;
	mat.NZ_size = NZ_size;
	mat.row_size = row_size;
	mat.col_size = column_size;
	mat.index = (coo_tuple*) malloc(NZ_size * sizeof(coo_tuple));
	return mat;
}


void coo_add_tuple(COO* matrix, long row, long col, long data, long index){
	coo_tuple new_tuple = { row, col, data };
	(*matrix).index[index] = new_tuple;
}


void reset_matrix(){
    long i,j;
    for(i=0;i<20;i++){
        for(j=0;j<20;j++)
            print_mat[i][j] = 0;
    }
}

void print_csr_matrix(CSR (*matrix)){
	long i,j, row_size, col_size, k, ind;
	row_size = (*matrix).row_size;
	col_size = (*matrix).col_size;
	//printf("row: %ld col: %ld\n", row_size, col_size);
	for(i=0;i<row_size;i++){
        j = (*matrix).row_ptr[i];
        k = (*matrix).row_ptr[i+1];
        for(ind=j-1;ind<k-1;ind++){
            print_mat[i][(*matrix).col_ind[ind]-1] = (*matrix).val[ind];
        }
	}
	for(i=0;i<row_size;i++){
		for(j=0;j<col_size;j++){
			printf("%ld\t", print_mat[i][j]);
		}
		printf("\n");
	}
	reset_matrix();
}

void print_csc_matrix(CSC (*matrix)){
	long i,j, row_size, col_size, k, ind;
	row_size = (*matrix).row_size;
	col_size = (*matrix).col_size;
	for(i=0;i<col_size;i++){
        j = (*matrix).col_ptr[i];
        k = (*matrix).col_ptr[i+1];
        for(ind=j-1;ind<k-1;ind++){
            print_mat[(*matrix).row_ind[ind]-1][i] = (*matrix).val[ind];
        }
	}
	for(i=0;i<row_size;i++){
		for(j=0;j<col_size;j++){
			printf("%ld\t", print_mat[i][j]);
		}
		printf("\n");
	}
	reset_matrix();
}

long find_node(COO_LL* head, long row, long col, long val){
    COO_LL *p=head;
    while(p != NULL){
        if((*p).coo_entry.row == row && (*p).coo_entry.col == col){
            (*p).coo_entry.val += val;
            return 1;
        }
        p = p->next;
    }
    return 0;
}
void printll(COO_LL* head){
    COO_LL *p=head;
    while(p != NULL){
        printf("row: %ld col: %ld val: %ld\t", p->coo_entry.row, p->coo_entry.col, p->coo_entry.val);
        p = p->next;
    }
    printf("\n");
    return ;
}
void add_node(COO_LL** ll, long row, long col, long val){
    //printf("row: %ld col: %ld val: %ld\n", row, col, val);
    if(*ll == NULL){
        COO_LL* new_node = (COO_LL*) malloc(sizeof(COO_LL));
        new_node->coo_entry.row = row;
        new_node->coo_entry.col = col;
        new_node->coo_entry.val = val;
        new_node->next = *ll;
        *ll = new_node;
        return ;
    }
    if(find_node(*ll, row, col, val)==0){
        COO_LL* new_node = (COO_LL*) malloc(sizeof(COO_LL));
        new_node->coo_entry.row = row;
        new_node->coo_entry.col = col;
        new_node->coo_entry.val = val;
        new_node->next = *ll;
        *ll = new_node;
    }
    /*printf("printing linked list after adding node\n");
    printll(*ll);*/
    return ;
}

COO create_coo_from_LL(COO_LL** ll, long row_size, long col_size){
    long NZ_size = 0, i = 0;
    COO_LL* p;
    p = *ll;
    while(p != NULL){
        if((*p).coo_entry.val != 0)
            NZ_size++;
        p = p->next;
    }
    COO coo_result = create_coo_matrix(row_size, col_size, NZ_size);
    p=*ll;
    while(p != NULL){
        if((*p).coo_entry.val != 0){
            coo_result.index[i++] = (*p).coo_entry;
        }
        p = p->next;
    }
    return coo_result;
}

COO csr_mul_csc(CSR* csr_mat, CSC* csc_mat){
	COO coo_mat;
    if(csr_mat->col_size != csc_mat->row_size){
        printf("Invalid Matrix Multiplication\n");
		coo_mat.NZ_size = -1;
        return coo_mat;
    }
    long i, j = 0, k, coo_len = 0, l;
    COO_LL* coo_ll = NULL;
    while(csr_mat->row_ptr[j] == -1){
        j++;
    }
    for(i=0;i<csr_mat->NZ_size;i++){
        if(csc_mat->col_ptr[j+1] == i+1){
            j++;
        }
        l=0;
        while(csc_mat->col_ptr[l] == -1){
            l++;
        }
        for(k=0;k<csc_mat->NZ_size;k++){
            if(csc_mat->col_ptr[l+1] == k+1){
                l++;
            }
            if(csr_mat->col_ind[i] == csc_mat->row_ind[k]){
                add_node(&coo_ll, j+1, l+1, csr_mat->val[i]*csc_mat->val[k]);
            }
        }
    }
    coo_mat = create_coo_from_LL(&coo_ll, csr_mat->row_size, csc_mat->col_size);
	return coo_mat;
}

COO csc_mul_csr(CSC* csc_mat, CSR* csr_mat){
	COO coo_mat;
    if(csc_mat->col_size != csr_mat->row_size){
        printf("Invalid Matrix Multiplication\n");
		coo_mat.NZ_size = -1;
        return coo_mat;
    }
    long i, j = 0, k, coo_len = 0, l=0;
    COO_LL* coo_ll = NULL;

    while(csc_mat->col_ptr[l] == -1){
            l++;
        }
    for(i=0;i<csc_mat->NZ_size;i++){

        if(csc_mat->col_ptr[l+1] == i+1){
            l++;
        }
        j=0;
        while(csr_mat->row_ptr[j] == -1){
            j++;
        }
        for(k=0;k<csr_mat->NZ_size;k++){
            if(csc_mat->col_ptr[j+1] == k+1){
                j++;
            }
            if(l == j){
                add_node(&coo_ll, csc_mat->row_ind[i], csr_mat->col_ind[k], csc_mat->val[i]*csr_mat->val[k]);
            }
        }
    }
    coo_mat = create_coo_from_LL(&coo_ll, csc_mat->row_size, csr_mat->col_size);
	return coo_mat;
}

void print_coo_matrix(COO (*matrix)){
	long i,j;
	for(i=0;i<(*matrix).NZ_size;i++){
		print_mat[(*matrix).index[i].row-1][(*matrix).index[i].col-1] = (*matrix).index[i].val;
	}
	for(i=0;i<(*matrix).row_size;i++){
		for(j=0;j<(*matrix).col_size;j++){
			printf("%ld\t", print_mat[i][j]);
		}
		printf("\n");
	}
	reset_matrix();
}

void swap(coo_tuple* a, coo_tuple* b){
    coo_tuple t = *a;
    *a = *b;
    *b = t;
}

long partition_row(coo_tuple coo_arr[], long low, long high){
    coo_tuple pivot = coo_arr[high];
    long i = (low - 1), j;
    for (j = low; j <= high- 1; j++)
    {
        if ((coo_arr[j].row < pivot.row) || ((coo_arr[j].row == pivot.row) && (coo_arr[j].col < pivot.col)))
        {
            i++;
            swap(&coo_arr[i], &coo_arr[j]);
        }
    }
    swap(&coo_arr[i + 1], &coo_arr[high]);
    return (i + 1);
}

long partition_col(coo_tuple coo_arr[], long low, long high){
    coo_tuple pivot = coo_arr[high];
    long i = (low - 1), j;
    for (j = low; j <= high- 1; j++)
    {
        if ((coo_arr[j].col < pivot.col) || ((coo_arr[j].col == pivot.col) && (coo_arr[j].row < pivot.row)))
        {
            i++;
            swap(&coo_arr[i], &coo_arr[j]);
        }
    }
    swap(&coo_arr[i + 1], &coo_arr[high]);
    return (i + 1);
}

void coo_quickSort_col(coo_tuple coo_arr[], long low, long high){
    if (low < high)
    {
        long part_index = partition_col(coo_arr, low, high);
        coo_quickSort_col(coo_arr, low, part_index - 1);
        coo_quickSort_col(coo_arr, part_index + 1, high);
    }
}

void coo_quickSort_row(coo_tuple coo_arr[], long low, long high){
    if (low < high)
    {
        long part_index = partition_row(coo_arr, low, high);
        coo_quickSort_row(coo_arr, low, part_index - 1);
        coo_quickSort_row(coo_arr, part_index + 1, high);
    }
}

void coo_sort_col(COO* coo_mat){
	coo_quickSort_col((*coo_mat).index, 0, (*coo_mat).NZ_size - 1);
}

CSC create_csc_from_coo(COO* coo_mat){
	CSC csc_mat = create_csc((*coo_mat).NZ_size, (*coo_mat).row_size, (*coo_mat).col_size);
	//printf("row: %ld col: %ld\n", csc_mat.row_size, csc_mat.col_size);
	long i,j=0, col_ind = 1, flag = 0;
	for(i=0;i<=(*coo_mat).col_size;i++){
		csc_mat.col_ptr[i] = -1;
	}
	coo_sort_col(coo_mat);

	coo_tuple temp;
	temp = (*coo_mat).index[0];
	csc_mat.val[0] = temp.val;
	csc_mat.row_ind[0] = temp.row;
	csc_mat.col_ptr[temp.col-1] = 1;

	for(i=0;i<(*coo_mat).NZ_size;i++){
		temp = (*coo_mat).index[i];
		csc_mat.val[i] = temp.val;
		csc_mat.row_ind[i] = temp.row;
		if(temp.col != col_ind){
			csc_mat.col_ptr[temp.col-1] = i+1;
			col_ind = temp.col;
		}
	}
	if(csc_mat.col_ptr[(*coo_mat).col_size] == -1)
		csc_mat.col_ptr[(*coo_mat).col_size] = (*coo_mat).NZ_size + 1;
	for(i=(*coo_mat).col_size;i>=0;i--){
		if(csc_mat.col_ptr[i] == 1)
			break;
		if(csc_mat.col_ptr[i] == -1)
			csc_mat.col_ptr[i] = csc_mat.col_ptr[i+1];
	}
    return csc_mat;
}

void coo_sort_row(COO* coo_mat){
	coo_quickSort_row((*coo_mat).index, 0, (*coo_mat).NZ_size - 1);
}

CSR create_csr_from_coo(COO* coo_mat){
	CSR csr_mat = create_csr((*coo_mat).NZ_size, (*coo_mat).row_size, (*coo_mat).col_size);
	//printf("row: %ld col: %ld\n", csr_mat.row_size, csr_mat.col_size);
	long i,j=0, row_ind = 1, flag = 0;
	for(i=0;i<=(*coo_mat).row_size;i++){
		csr_mat.row_ptr[i] = -1;
	}
	coo_sort_row(&(*coo_mat));

	coo_tuple temp;
	temp = (*coo_mat).index[0];
	csr_mat.val[0] = temp.val;
	csr_mat.col_ind[0] = temp.col;
	csr_mat.row_ptr[temp.row-1] = 1;

	for(i=0;i<(*coo_mat).NZ_size;i++){
		temp = (*coo_mat).index[i];
		csr_mat.val[i] = temp.val;
		csr_mat.col_ind[i] = temp.col;
		if(temp.row != row_ind){
			csr_mat.row_ptr[temp.row-1] = i+1;
			row_ind = temp.row;
		}
	}
	if(csr_mat.row_ptr[(*coo_mat).row_size] == -1)
		csr_mat.row_ptr[(*coo_mat).row_size] = (*coo_mat).NZ_size + 1;
	for(i=(*coo_mat).row_size;i>=0;i--){
		if(csr_mat.row_ptr[i] == 1)
			break;
		if(csr_mat.row_ptr[i] == -1)
			csr_mat.row_ptr[i] = csr_mat.row_ptr[i+1];
	}
    return csr_mat;
}

void csc_row_transform(CSC* csc_mat, long r1, long k1, long r2, long k2){
    long i, j = 0, row1_arr[(*csc_mat).col_size], row2_arr[(*csc_mat).col_size], new_size = 0;
    for(i=0;i<(*csc_mat).col_size;i++){
        row1_arr[i] = 0;
        row2_arr[i] = 0;
    }
    while(csc_mat->col_ptr[j] == -1){
        j++;
    }
    for(i=0;i<csc_mat->NZ_size;i++){
        if(csc_mat->col_ptr[j+1] == i+1){
            j++;
        }
        if((*csc_mat).row_ind[i] == r1){
            row1_arr[j] = csc_mat->val[i];
        }
        if((*csc_mat).row_ind[i] == r2){
            row2_arr[j] = csc_mat->val[i];
        }
        else
            new_size++;
    }
    for(i=0;i<(*csc_mat).col_size;i++){
        row2_arr[i] = k2*row2_arr[i] + k1*row1_arr[i];
        if(row2_arr[i] != 0)
            new_size++;
    }
    COO coo_mat= create_coo_matrix(csc_mat->row_size, csc_mat->col_size, new_size);
    j=0;
    long ind=0;
    for(i=0;i<csc_mat->NZ_size;i++){
        if(csc_mat->col_ptr[j+1] == i+1){
            j++;
        }
        if((*csc_mat).row_ind[i] != r2){
            coo_tuple new_tuple = {(*csc_mat).row_ind[i], j+1, (*csc_mat).val[i]};
            coo_mat.index[ind++] = new_tuple;
        }
    }
    for(i=0;i<(*csc_mat).col_size;i++){
        if(row2_arr[i] != 0){
            coo_tuple new_tuple = {r2, i+1, row2_arr[i]};
            coo_mat.index[ind++] = new_tuple;
        }
    }
    *csc_mat = create_csc_from_coo(&coo_mat);
}

long check_cell_indices(COO* coo_mat, long row_ind, long col_ind, long size){
    long i;
    //printf("row: %ld col: %ld\n", row_ind, col_ind);
    for(i=0;i<size;i++){
        if((*coo_mat).index[i].row == row_ind && (*coo_mat).index[i].col == col_ind)
            return 0;
    }
    return 1;
}

COO create_random_sparse_matrix(long m, long n, long elems, long neg_nums){
    long i, j=0, row, column, data, MOD = 100;
    //printf("elems: %ld\n", elems);
	COO newMat = create_coo_matrix(m, n, elems);
	srand(time(0));
	for(i=0;i<elems;i++){
        row = rand() % m;
        column = rand() % n;
        //printf("row: %ld col: %ld\n", row, column);
        while(1){
            if(check_cell_indices(&newMat, row+1, column+1, i) == 1)
                break;
            row = rand() % m;
            column = rand() % n;
        }
        data = rand() % MOD;
        if(j<neg_nums){
            data *= -1;
            j++;
        }
		coo_add_tuple(&newMat, row + 1, column + 1, data + 1, i);
	}
	return newMat;
}
