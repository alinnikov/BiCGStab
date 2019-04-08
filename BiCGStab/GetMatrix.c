#include "GetMatrix.h"



int get_num_rows(char *name) {

	FILE *matrix_file;
	char str[1024];
	int count_of_values = 0;
	int num_rows = 0;

	if ((matrix_file = fopen(name, "rb")) == NULL) {
		printf("Cannot open file.\n");
	}
	//���������� ������ 2 ������
	for (int i = 0; i < 2; i++)
	{
		fgets(str, 1024, matrix_file);
		//printf("%s", str);
	}
	//���� ����� �����
	while ((fscanf(matrix_file, "%s", &str) != EOF && count_of_values < 3))
	{
		if (!matrix_file) break;    //����� �� ����� �������
		count_of_values++;
		if (count_of_values == 2) sscanf(str, "%d", &num_rows);
	}
	fclose(matrix_file);
	return num_rows;
}

int get_num_values(char *name) {
	FILE *matrix_file;
	char str[1024];
	int count_of_values = 0;
	int num_rows = 0;

	if ((matrix_file = fopen(name, "rb")) == NULL) {
		printf("Cannot open file.\n");
	}
	//���������� ������ 2 ������
	for (int i = 0; i < 2; i++)
	{
		fgets(str, 1024, matrix_file);
		//printf("%s", str);
	}
	//���� ����� ��������
	while ((fscanf(matrix_file, "%s", &str) != EOF && count_of_values < 5))
	{
		if (!matrix_file) break;    //����� �� ����� �������
		count_of_values++;
		if (count_of_values == 4) sscanf(str, "%d", &num_rows);
	}
	fclose(matrix_file);
	return num_rows;
}

//�������� CSR-�������
void get_csr_matrix(struct CSR_matrix *m, char *name) {

	//������� ���������� ���������
	m->num_rows = get_num_rows(name);
	m->num_values = get_num_values(name);
	m->array_rows = (int*)malloc((m->num_rows + 1) * sizeof(int));
	m->array_columns = (int*)malloc((m->num_values) * sizeof(int));
	m->array_values = (double*)malloc((m->num_values + 1) * sizeof(double));

	//��������� ���� � ��������
	FILE *fp;
	char str[1024];
	int k = 0;
	if ((fp = fopen(name, "rb")) == NULL) {
		printf("Cannot open file.\n");
	}
	//���������� ������ � ����������� ����������� � �������
	for (int i = 0; i < 4; i++)
	{
		fgets(str, 1024, fp);
	}

	//�������� ������ �����
	while (!feof(fp) && k < m->num_rows + 1) {


		if (fscanf(fp, "%s", str))
			sscanf(str, "%d", &m->array_rows[k]);
		k++;
	}

	//�������� ������ ��������
	k = 0;
	while (!feof(fp) && k < m->num_values) {
		if (fscanf(fp, "%s", str))
			sscanf(str, "%d", &m->array_columns[k]);
		//printf("\n%d", array2[k]);
		k++;
	}

	//�������� ������ ��������
	k = 0;
	while (!feof(fp)) {
		if (fscanf(fp, "%s", str))
			sscanf(str, "%lf", &m->array_values[k]);
		//printf("\n%f", array3[k]);
		k++;
	}

	//printf("%lf", m->array_values[k - 1]);
}