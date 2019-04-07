#include "GetMatrix.h"



int get_num_rows() {

	FILE *matrix_file;
	char str[1024];
	int count_of_values = 0;
	int num_rows1 = 0;

	if ((matrix_file = fopen("matrix.rb", "rb")) == NULL) {
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
		if (count_of_values == 2) sscanf(str, "%d", &num_rows1);
	}
	fclose(matrix_file);
	return num_rows1;
}

int get_num_values() {
	FILE *matrix_file;
	char str[1024];
	int count_of_values = 0;
	int num_rows = 0;

	if ((matrix_file = fopen("matrix.rb", "rb")) == NULL) {
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
void get_csr_matrix(struct CSR_matrix *m) {

	//������� ���������� ���������
	m->num_rows = get_num_rows();
	m->num_values = get_num_values();
	m->array_rows = (int*)malloc((m->num_rows + 1) * sizeof(int));
	m->array_columns = (int*)malloc((m->num_values) * sizeof(int));
	m->array3 = (float*)malloc((m->num_values + 1) * sizeof(float));

	//��������� ���� � ��������
	FILE *fp;
	char str[1024];
	int k = 0;
	if ((fp = fopen("matrix.rb", "rb")) == NULL) {
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
			sscanf(str, "%f", &m->array3[k]);
		//printf("\n%f", array3[k]);
		k++;
	}

	printf("%f", m->array3[k - 1]);
}