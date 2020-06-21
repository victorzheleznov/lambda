#include <conio.h>
#include <stdio.h>
#include <stdlib.h>

#include "lambda.h"

int main(void)
{
	//double matr[9] = {6.290, 5.978, 0.544, 5.978, 6.292, 2.340, 0.544, 2.340, 6.288};
	//double trian_matr[9] = {0.0};
	//double diag_matr[3] = {0.0};
	//double est[12] = {-28490.8566886116, 65752.6299198198, 38830.3666554972, 5003.70833517778, -29196.0699104593, -297.658932458787, -22201.0284440701, 51235.8374755528, 30257.7809603224, 3899.40332138829, -22749.1853575113, -159.278779870217};

	candidates* cands = {0};

	int dim = 6;
	int num_cands = 2;
	int found = 0;

	int* temp = NULL;

	double trian_matr[36] = {1, 0, 0, 0, 0, 0, 
							 0.20133153553671487, 1, 0, 0, 0, 0, 
		                    -0.40263455477678778, 0.22836520497039933, 1, 0, 0, 0, 
		                     0.53996647368506812, 0.42120362186266830, -0.059202661566631468, 1, 0, 0, 
		                     0.095522480284414249, 0.87426776631741099, 0.89110797471714276, 0.62618725559123201, 1, 0, 
		                     0.95842519013860772, -0.51375074226063455, -1.5270588217087484, 0.068501493402481073, -0.73074292414498698, 1};
	
	
	double diag_matr[6] = {0.00083662115115058658, 0.00072600689290804737, 0.0014435334490503994, 0.0010231785948171545, 0.0035376786776778820, 0.0015875936529177563};
	
	double est[6] = {-101.92854683847105, 71.000646234118349, -5.0909824324773982, 111.02821414315092, 117.97630226915967, -85.949298535284598};
	int transf_matr[36] = {0};

	FILE* transformation_matrix_txt = NULL;
	FILE* candidates_txt = NULL;
		
	fopen_s(&transformation_matrix_txt, "transformation_matrix.txt", "w");
	fopen_s(&candidates_txt, "candidates.txt", "w");

	transformation(dim, trian_matr, diag_matr, est, transf_matr);

	for (int i = 0; i < dim; i++)
	{
		for (int j = 0; j < dim; j++)
		{
			fprintf(transformation_matrix_txt, "%-5d", transf_matr[i * dim + j]);
		}
		fprintf(transformation_matrix_txt, "\n");
	}


	cands = (candidates*) malloc(num_cands * sizeof(candidates));
	
	for (int j = 0; j < num_cands; j++)
	{
		cands[j].est = NULL;
		cands[j].norm = 0.0;

		cands[j].est = (int*) malloc(dim * sizeof(int));
	}

	search(dim, num_cands, trian_matr, diag_matr, est, &found, cands);

	// back transformation

	temp = (int*) malloc(dim * sizeof(int));

	for (int k = 0; k < found; k++)
	{
		for (int i = 0; i < dim; i++)
		{
			temp[i] = 0;
			for (int j = 0; j < dim; j++) temp[i] += transf_matr[i * dim + j] * cands[k].est[j];
		}

		for (int i = 0; i < dim; i++) cands[k].est[i] = temp[i];
	}

	// output

	for (int i = 0; i < num_cands; i++)
	{
		for (int j = 0; j < dim; j++)
		{
			fprintf(candidates_txt, "%-5d", cands[i].est[j]);
		}
		fprintf(candidates_txt, "\n");
	}

	//_getch();

	for (int i = 0; i < num_cands; i++)
	{
		free(cands[i].est);
	}

	free(cands);

	fclose(transformation_matrix_txt);
	fclose(candidates_txt);
}