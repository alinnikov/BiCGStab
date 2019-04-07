#include "GetMatrix.h"



int get_num_rows() {

	FILE *matrix_file;
	char str[1024];
	int count_of_values = 0;
	int num_rows = 0;

	if ((matrix_file = fopen("matrix.rb", "rb")) == NULL) {
		printf("Cannot open file.\n");
	}
	//Пропускаем первые 2 строки
	for (int i = 0; i < 2; i++)
	{
		fgets(str, 1024, matrix_file);
		//printf("%s", str);
	}
	//Ищем число строк
	while ((fscanf(matrix_file, "%s", &str) != EOF && count_of_values < 3))
	{
		if (!matrix_file) break;    //чтобы не делал лишнего
		count_of_values++;
		if (count_of_values == 2) sscanf(str, "%d", &num_rows);
	}
	fclose(matrix_file);
	return num_rows;
}

int get_num_values() {
	FILE *matrix_file;
	char str[1024];
	int count_of_values = 0;
	int num_rows = 0;

	if ((matrix_file = fopen("matrix.rb", "rb")) == NULL) {
		printf("Cannot open file.\n");
	}
	//Пропускаем первые 2 строки
	for (int i = 0; i < 2; i++)
	{
		fgets(str, 1024, matrix_file);
		//printf("%s", str);
	}
	//Ищем число значений
	while ((fscanf(matrix_file, "%s", &str) != EOF && count_of_values < 5))
	{
		if (!matrix_file) break;    //чтобы не делал лишнего
		count_of_values++;
		if (count_of_values == 4) sscanf(str, "%d", &num_rows);
	}
	fclose(matrix_file);
	return num_rows;
}

//Получаем CSR-матрицу
void get_csr_matrix(struct CSR_matrix *m) {

	//Создаем переменные структуры
	m->num_rows = get_num_rows();
	m->num_values = get_num_values();
	m->array_rows = (int*)malloc((m->num_rows + 1) * sizeof(int));
	m->array_columns = (int*)malloc((m->num_values) * sizeof(int));
	m->array_values = (float*)malloc((m->num_values + 1) * sizeof(float));

	//Открываем файл с матрицей
	FILE *fp;
	char str[1024];
	int k = 0;
	if ((fp = fopen("matrix.rb", "rb")) == NULL) {
		printf("Cannot open file.\n");
	}

	//Пропускаем строки с технической информацией о матрице
	for (int i = 0; i < 4; i++)
	{
		fgets(str, 1024, fp);
	}

	//Получаем массив строк
	while (!feof(fp) && k < m->num_rows + 1) {


		if (fscanf(fp, "%s", str))
			sscanf(str, "%d", &m->array_rows[k]);
		k++;
	}

	//Получаем массив столбцов
	k = 0;
	while (!feof(fp) && k < m->num_values) {
		if (fscanf(fp, "%s", str))
			sscanf(str, "%d", &m->array_columns[k]);
		//printf("\n%d", array2[k]);
		k++;
	}

	//Получаем массив значений
	k = 0;
	while (!feof(fp)) {
		if (fscanf(fp, "%s", str))
			sscanf(str, "%f", &m->array_values[k]);
		//printf("\n%f", array3[k]);
		k++;
	}

	//printf("%f", m->array3[k - 1]);
}