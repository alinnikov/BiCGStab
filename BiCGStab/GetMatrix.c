#include "GetMatrix.h"


int get_num_rows(char *name) {

	FILE *matrix_file;
	char str[1024];
	int count_of_values = 0;
	int num_rows = 0;

	if ((matrix_file = fopen(name, "rb")) == NULL) {
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
	return num_rows-1;
}

int get_num_values(char *name) {
	FILE *matrix_file;
	char str[1024];
	int count_of_values = 0;
	int num_values = 0;

	if ((matrix_file = fopen(name, "rb")) == NULL) {
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
		if (count_of_values == 4) sscanf(str, "%d", &num_values);
	}
	fclose(matrix_file);
	return num_values-1;
}

//Получаем CSR-матрицу
void get_csr_matrix(struct CSR_matrix *m, char *name) {

	//Создаем переменные структуры
	m->num_rows = get_num_rows(name);
	m->num_values = get_num_values(name);
	m->array_rows = (int*)malloc((m->num_rows+1) * sizeof(int));
	m->array_columns = (int*)malloc((m->num_values) * sizeof(int));
	m->array_values = (double*)malloc((m->num_values) * sizeof(double));
	//Открываем файл с матрицей
	FILE *fp;
	char str[1024];
	int k = 0;
	if ((fp = fopen(name, "rb")) == NULL) {
		printf("Cannot open file.\n");
	}
	//Пропускаем строки с технической информацией о матрице
	for (int i = 0; i < 4; i++)
	{
		fgets(str, 1024, fp);
	}

	//Получаем массив строк
	while (!feof(fp) && k <= m->num_rows+1) {
		if (fscanf(fp, "%s", str)) {
			sscanf(str, "%d", &m->array_rows[k]);
			m->array_rows[k]--;
		}
		k++;
	}
	//Получаем массив столбцов
	k = 0;
	while (!feof(fp) && k <= m->num_values) {
		if (fscanf(fp, "%s", str)) {
			sscanf(str, "%d", &m->array_columns[k]);
			m->array_columns[k]--;
		}
		k++;
	}

	//Получаем массив значений
	k = 0;
	while (!feof(fp)) {
		if (fscanf(fp, "%s", str)) {
			sscanf(str, "%lf", &m->array_values[k]);
		}
		k++;
	}
}