#ifndef ILU_H
#define ILU_H

#include "linalg.h"
#include <vector>

/*
 *  Perform incomplete LU factorization of A, where A is stored in
 *  CSR format.
 *
 *  Parameters
 *      Ap       - CSR row poIer
 *      Aj       - CSR index array
 *      Ax       - CSR data array
 * 		Sp       - CSR row pointer
 *      Sj       - CSR index array
 *      Sx       - CSR data array
 *
 *  Returns:
 *      Nothing, Sp, Sj and Sx will be modified in place
 *
 */
template<class I, class T>
void ilu_levels(const I Ap[], const int Ap_size,
                  const I Aj[], const int Aj_size, 
                  const T Ax[], const int Ax_size, 
                  I Sp[], int Sp_size,
                  I Sj[], int Sj_size,
                  T Sx[], int Sx_size)
{
	I n_rows = Ap_size-1;
	I levls_nnz = 0;

	Sp[0] = 0;

	I start_r = Ap[0];
 	I end_r = Ap[1];
 	for(I i = start_r; i<end_r;i++){
          Sj[levls_nnz] = Aj[i];
          Sx[levls_nnz] = 0;
          levls_nnz++;
    }
    Sp[1] = levls_nnz;

    std::vector<I> current_row_levls;
    std::vector<I> current_row_p;
    std::vector<I> temp_row_levls;
    std::vector<I> temp_row_p;

    current_row_levls.reserve(n_rows);
    current_row_p.reserve(n_rows);
    temp_row_levls.reserve(n_rows);
    temp_row_p.reserve(n_rows);

	for(I row_i = 1; row_i < n_rows;row_i++){
		I start_ri = Ap[row_i];
		I end_ri = Ap[row_i+1];
		I nnz_i = 0;
		I temp_nnz_i = 0;
		I current_row_ind = 0;

		//initialize temporary levls and idx2 vectors for row i
		for(I j = start_ri; j < end_ri; j++){
			current_row_levls.push_back(0);
			current_row_idx2.push_back(Aj[j]);
			nnz_i++;
		}

		for(I j = start_ri; j < end_ri; j++){
			I col_k = Aj[j];
			
			if(col_k >= row_i){
				break;
			}
			
			Sx[levls_nnz] = (current_row_levls[current_row_ind]);
			Sj[levls_nnz] = (current_row_idx2[current_row_ind]);
			
			current_row_ind++;
			levls_nnz++;		


			//get row k from updated levls matrix
			I start_rk = Sp[col_k];
			I end_rk = Sp[col_k+1];
			
			//temp variables
			I ind_i = -1;
			I ind_k = -1;
			
			I col_ij = -1;
			I col_kj = -1;
			I col_ik = col_k;

			I lev_ij = 100000;
			I lev_ik = 100000;
			I lev_kj = 100000;
			
			//get starting indices for rows i and k such that col >=k+1
			for(I i = 0; i < nnz_i; i++){
				I col = current_row_idx2[i];
				if(col == col_k)
					lev_ik = current_row_levls[i];
				if(col >= col_k+1){
					ind_i = i;
					col_ij = current_row_idx2[i];
					break;
				}
			}


			for(I i = start_rk; i < end_rk; i++){
				I col = Sp[i];
				if(col >=col_k+1){
					ind_k = i;
					col_kj = Sj[i];
					break;
				}
			}

			//check if reached end of row k
			if(ind_k == -1)
				col_kj = n_rows + 1;

			I current_col = -1;
			I current_row = -1;

			while(1){
				if((col_ij >= n_rows) and (col_kj >= n_rows)){
					break;
				}

				if(col_ij == col_kj){
					current_col = col_ij;
					current_row = row_i;

					lev_ij = current_row_levls[ind_i];
					lev_kj = Sx[ind_k];

					lev_ij = min(lev_ij, lev_ik+lev_kj+1);

					if(lev_ij < 10000){
						temp_row_levls.push_back(lev_ij);
						temp_row_idx2.push_back(current_col);
						temp_nnz_i++;
					}

					ind_i++;
					ind_k++;
				}

				else if(col_ij<col_kj){
					current_col=col_ij;
					current_row=row_i;

					lev_ij = current_row_levls[ind_i];

					if(lev_ij < 10000){
						temp_row_levls.push_back(lev_ij);
						temp_row_idx2.push_back(current_col);
						temp_nnz_i++;
					}

					ind_i++;
				}

				else if(col_ij>col_kj){
					current_col = col_kj;
					current_row = col_k;
					
					lev_kj = Sx[ind_k];

					lev_ij = lev_ik + lev_kj + 1;
					
					if(lev_ij < 10000){
						temp_row_levls.push_back(lev_ij);
						temp_row_idx2.push_back(current_col);
						temp_nnz_i++;
					}

					ind_k++;
				}

				if((ind_i < nnz_i) && (ind_i != -1))
					col_ij=current_row_idx2[ind_i];
				else
					col_ij = n_rows+1;

				if((ind_k < end_rk) && (ind_k != -1))
					col_kj=Sj[ind_k];
				else
					col_kj = n_rows+1;
			}
			
			current_row_levls = temp_row_levls;
			current_row_idx2 = temp_row_idx2;
			nnz_i = temp_nnz_i;
			
			temp_row_levls.clear();
			temp_row_idx2.clear();
			temp_nnz_i = 0; 	
		}
		
		//Add current row levels to final levls matrix
		for(I i = 0; i<nnz_i;i++){
			Sj[levls_nnz] = (current_row_idx2[i]);
			Sx[levls_nnz] = (current_row_levls[i]);
			levls_nnz++;
		}

		current_row_levls.clear();
		current_row_idx2.clear();

		Sp[row_i+1] = levls_nnz;
	}//end i loop
}

/*
 *  Gets desired sparsity pattern of ILU factors based 
 *  on level of fill.
 *
 *  Parameters
 *      Ap       - CSR row pointer of levels of fill matrix
 *      Aj       - CSR index array of levels of fill matrix
 *      Ax       - CSR data array of levels of fill matrix
 * 		Sp       - CSR row pointer of sparsity matrix
 *      Sj       - CSR index array of sparsity matrix
 *      Sx       - CSR data array of sparsity matrix
 *      lof      - Level of fill
 *
 *  Returns:
 *      Nothing, Sp, Sj and Sx will be modified in place
 *
 */
template<class I, class T>
void ilu_sparsity(const I Ap[], const int Ap_size,
                  const I Aj[], const int Aj_size, 
                  const T Ax[], const int Ax_size, 
                  I Sp[], int Sp_size,
                  I Sj[], int Sj_size,
                  T Sx[], int Sx_size,
                  int lof)
{

	I Sn_rows = Ap_size-1;
	I Snnz = 0;
	
	Sp[0] = 0;
	
	for(I row_i = 0; row_i<Sn_rows;row_i++){
		Sp[row_i] = Snnz;
		I start_ri = Ap[row_i];
		I end_ri = Ap[row_i+1];
		for(I j = start_ri; j<end_ri;j++){
			I col_j = Aj[j];
			I levl_j = Ax[j];
			if(levl_j <= lof){
				Sj[Snnz] = (col_j);
				Sx[Snnz] = (levl_j);
				Snnz++;
			}
		}
	}
	Sp[Sn_rows] = Snnz;
}

/*
 *  Perform numeric phase of ILU factorization of A, where A is stored in
 *  CSR format.
 *
 *  Parameters
 *      Ap       - CSR row pointer
 *      Aj       - CSR index array
 *      Ax       - CSR data array
 * 		Sp       - CSR row pointer of sparsity
 *      Sj       - CSR index array of sparsity
 *      Sx       - CSR data array of sparsity
 *      Tx       - CSR data array of factors
 *
 *  Returns:
 *      Nothing, Tp, Tj and Tx will be modified in place
 *
 */
template<class I, class T>
void ilu_numeric(const I Ap[], const int Ap_size,
                  const I Aj[], const int Aj_size, 
                  const T Ax[], const int Ax_size, 
                  I Sp[], int Sp_size,
                  I Sj[], int Sj_size,
                  T Sx[], int Sx_size,
                  T Tx[], int Tx_size)
{
	I n_rows = Ap_size - 1;
	
	I index_i = 0;
	for(I i = 0; i<Sx_size;i++){
		if(Sx[i] == 0.0){
			Tx[i] = Ax[index_i];
			index_i++;
		}
		else
			Tx[i] = 0.0;
	}

	//get vector of diagonal elements
	std::vector<T> diag_vec(n_rows,0.0);
	
	for(I row = 0; row < n_rows; row++){
		I start_i = Sp[row];
		I end_i = Sp[row+1];
		for(I j = start_i; j<end_i; j++){
			I col = Sj[j];
			T val = Tx[j];
			if(row == col)
				diag_vec[row] = val;
		}
	}
	
	//ikj Gaussian Elimination
	for(I row_i = 1; row_i<n_rows; row_i++){
		I start_ri = Sp[row_i];
		I end_ri = Sp[row_i+1];

		for(I j = start_ri; j<end_ri; j++){
			I col_k = Sj[j];

			//since k loop goes from 0 to i-1
			if(col_k >= row_i){
				break;
			}

			//compute multiplier
			Tx[j] /= diag_vec[col_k];
			T mult = Tx[j];

			if(row_i == col_k){
				diag_vec[row_i] = mult;
			}
		
			//get row k from original matrix
			I start_rk = Ap[col_k];
			I end_rk = Ap[col_k+1];

			//temp variables
			I ind_i = -1;
			I ind_k = -1;

			I col_ij = -1;
			I col_kj = -1;
			I col_ik = col_k;
			
			for(I i=start_ri; i< end_ri;i++){
				I col = Sj[i];
				if(col >= col_k+1){
					ind_i = i;
					col_ij = Sj[i];
					break;
				}		
			}

			for(I i=start_rk; i<end_rk;i++){
				I col = idx2[i];
				if(col >= col_k+1){
					ind_k = i;
					col_kj = Sj[i];
					break;
				}
			}

			//check if reached end of row k
			if(ind_k == -1)
				col_kj = n_rows + 1;

			T current_it = 0.0;
			I current_col = -1;
			I current_row = -1;

			while(1){
				if((col_ij >= n_rows) && (col_kj >= n_rows))
					break;

				if(col_ij == col_kj){
					current_col = col_ij;
					current_row = row_i;
					current_it = Tx[ind_i] - mult*Tx[ind_k];
					Tx[ind_i] = current_it;
				
					//check if it's a diagonal element
					if(current_row == current_col)
						diag_vec[current_row] = current_it;
				
					ind_i++;
					ind_k++;	
				
				}	
				else if(col_ij<col_kj){
					current_col = col_ij;
					current_row = row_i;
					current_it = Tx[ind_i];
					Tx[ind_i] = current_it;

					//check if it's a diagonal element
					if(current_row == current_col)
						diag_vec[current_row] = current_it;

					ind_i++;
				}
				else if(col_ij>col_kj){
					current_col = col_kj; 
					current_row = col_k;
					current_it = -mult*Tx[ind_k];
					Tx[ind_k] = current_it;

					//check if it's a diagonal element
					if(current_row == current_col)
						diag_vec[current_row] = current_it;

					ind_k++;
				}

				if((ind_i<end_ri)&&(ind_i !=-1))
					col_ij = Sj[ind_i];
				else
					col_ij = n_rows + 1;

				if((ind_k<end_rk)&&(ind_k != -1))
					col_kj = Sj[ind_k];
				else
					col_kj = n_rows + 1;
			}
		}
	}
}






