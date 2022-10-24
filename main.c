#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct data_structure {
    unsigned magic_number;
    unsigned data_count;
    unsigned height, length;
    char **data;
} data_structure;

typedef struct label_structure {
    unsigned magic_number;
    unsigned data_count;
    char **data;
} label_structure;

unsigned reverse_int(unsigned rint) {
    return (unsigned) (((rint << 24) & 0xff000000) | ((rint << 8) & 0x00ff0000) |
                       ((rint >> 8) & 0x0000ff00) | ((rint >> 24) & 0x000000ff));
}

data_structure *load_data(const char *file_name) {
    FILE *data_file = fopen(file_name, "rb");
    data_structure *new_data_structure = (data_structure *) malloc(sizeof(data_structure));
    fread(new_data_structure, 4, 4, data_file);
    new_data_structure->data_count = reverse_int(new_data_structure->data_count);
    new_data_structure->height = reverse_int(new_data_structure->height);
    new_data_structure->length = reverse_int(new_data_structure->length);
    printf("load data file(%s) : data count = %u height = %u length = %u dimension = %d\n",
           file_name, new_data_structure->data_count, new_data_structure->height, new_data_structure->length,
           new_data_structure->height * new_data_structure->length);
    new_data_structure->data = (char **) malloc(sizeof(char *) * new_data_structure->data_count);
    unsigned data_depth = new_data_structure->height * new_data_structure->length;
    for (int i = 0; i < new_data_structure->data_count; ++i) {
        new_data_structure->data[i] = (char *) malloc(data_depth);
        fread(new_data_structure->data[i], 1, data_depth, data_file);
    }
    fclose(data_file);
    return new_data_structure;
}

label_structure *load_label(const char *file_name) {
    FILE *label_file = fopen(file_name, "rb");
    label_structure *new_label_structure = (label_structure *) malloc(sizeof(label_structure));
    fread(new_label_structure, 4, 2, label_file);
    new_label_structure->data_count = reverse_int(new_label_structure->data_count);
    printf("load label file(%s) : data count = %u\n", file_name, new_label_structure->data_count);
    new_label_structure->data = (char **) malloc(sizeof(char *) * new_label_structure->data_count);
    char *tmp = (char *) malloc(new_label_structure->data_count);
    fread(tmp, 1, new_label_structure->data_count, label_file);

    // weight
    unsigned char thres[10] = {30, 1, 32, 25, 25, 34, 30, 26, 28, 27};

    for (int i = 0; i < new_label_structure->data_count; ++i) {
        new_label_structure->data[i] = (char *) calloc(10, 1), new_label_structure->data[i][tmp[i]] = thres[tmp[i]];
    }
    free(tmp);
    fclose(label_file);
    return new_label_structure;
}

double determinant(double **matrix, unsigned length) {
    printf("reduce to length = %d", length);
    if (length == 1) return matrix[0][0];
    double **calculate_matrix = (double **) malloc(sizeof(double *) * (length - 1));
    for (int i = 0; i < length - 1; ++i) calculate_matrix[i] = (double *) malloc(sizeof(double) * (length - 1));
    for (int i = 1; i < length; ++i) for (int j = 1; j < length; ++j) calculate_matrix[i - 1][j - 1] = matrix[i][j];
    double sum = matrix[0][0] * determinant(calculate_matrix, length - 1), inverse = -1;
    for (int i = 0; i < length - 1; ++i) {
        for (int j = 1; j < length; ++j) {
            calculate_matrix[j - 1][i] = matrix[j][i];
        }
        sum += matrix[0][i + 1] * determinant(calculate_matrix, length - 1) * inverse;
        inverse *= -1;
    }

    for (int i = 0; i < length - 1; ++i) free(calculate_matrix[i]);
    free(calculate_matrix);
    return sum;
}

void process_row_add(unsigned from_line, unsigned to_line, double factor,
                     double **matrix, double **reverse, unsigned length) {
    for (int i = 0; i < length; ++i) {
        matrix[to_line][i] += matrix[from_line][i] * factor;
        reverse[to_line][i] += reverse[from_line][i] * factor;
    }
}

void process_divide_line(unsigned line, double factor,
                         double **matrix, double **reverse, unsigned length) {
    for (int i = 0; i < length; ++i) {
        matrix[line][i] /= factor;
        reverse[line][i] /= factor;
    }
}

void reverse_matrix(double **matrix, double **reverse, unsigned length) {
    for (int i = 0; i < length; ++i) memset(reverse[i], 0, 8 * length), reverse[i][i] = 1;
    // front eliminate
    for (int column = 0; column < length; ++column) {
        for (int i = column; i < length; ++i) {
            for (int j = i + 1; j < length; ++j) {
                if (matrix[i][column] == 0) matrix[i][column] = 0.0000001;
                process_row_add(i, j, -matrix[j][column] / matrix[i][column], matrix, reverse, length);
            }
        }
        printf("subprocess : gaussian_elimination : process [ %d / %d ]\n", column + 1, length);
    }
    // back substitute
    for (int column = length - 1; column >= 0; --column) {
        for (int i = column; i >= 0; --i) {
            for (int j = i - 1; j >= 0; --j) {
                if (matrix[i][column] == 0) matrix[i][column] = 0.0000001;
                process_row_add(i, j, -matrix[j][column] / matrix[i][column], matrix, reverse, length);
            }
        }
        printf("subprocess : back_substitution : process [ %d / %d ]\n", length - column, length);
    }
    // deal with diagonal
    printf("subprocess : extract diagonal matrix.\n");
    for (int i = 0; i < length; ++i) process_divide_line(i, matrix[i][i], matrix, reverse, length);
}

double **calculate_network(data_structure *train_data, label_structure *label_data) {
    unsigned data_depth = train_data->height * train_data->length;

    // multiple AT A => AT*A
    printf("main process : starting multiple ( A_transpose, A ).\n");
    double **result_matrix_ATA = (double **) malloc(sizeof(double *) * data_depth);
    for (int i = 0; i < data_depth; ++i) result_matrix_ATA[i] = (double *) malloc(sizeof(double) * data_depth);
    for (int i = 0; i < data_depth; ++i) {
        for (int j = 0; j < data_depth; ++j) {
            unsigned sum = 0;
            for (int k = 0; k < train_data->data_count; ++k) {
                sum += ((unsigned) (unsigned char) train_data->data[k][i]) *
                       ((unsigned) (unsigned char) train_data->data[k][j]);
            }
            result_matrix_ATA[i][j] = (double) sum;
        }
        printf("subprocess : process [ %d / %d ]\n", i + 1, data_depth);
    }
    printf("main process : ended multiple ( A_transpose, A ).\n");

    // reverse (AT*A) => (AT*A)-1
    printf("main process : starting reverse ( A_transpose * A ).\n");
    double **result_cofactor_transpose = (double **) malloc(sizeof(double *) * data_depth);
    for (int i = 0; i < data_depth; ++i) result_cofactor_transpose[i] = (double *) malloc(sizeof(double) * data_depth);
    reverse_matrix(result_matrix_ATA, result_cofactor_transpose, data_depth);
    for (int i = 0; i < data_depth; ++i) free(result_matrix_ATA[i]);
    free(result_matrix_ATA);
    printf("main process : ended reverse ( A_transpose * A ).\n");

    // multiple AT b => AT*b
    printf("main process : starting multiple ( A_transpose, b ).\n");
    double **result_matrix_ATB = (double **) malloc(sizeof(double *) * data_depth);
    for (int i = 0; i < data_depth; ++i) result_matrix_ATB[i] = (double *) malloc(sizeof(double) * 10);
    for (int i = 0; i < data_depth; ++i) {
        for (int j = 0; j < 10; ++j) {
            unsigned sum = 0;
            for (int k = 0; k < train_data->data_count; ++k) {
                sum += ((unsigned) (unsigned char) train_data->data[k][i]) * ((unsigned) label_data->data[k][j]);
            }
            result_matrix_ATB[i][j] = sum;
        }
        printf("subprocess : process [ %d / %d ] \n", i + 1, data_depth);
    }
    printf("main process : ended multiple ( A_transpose, b ).\n");

    // multiple (AT*A)-1 AT*b => (AT*A)-1ATb => x
    printf("main process : starting multiple ( reverse ( A_transpose * A ), A * b ).\n");
    double **network = (double **) malloc(sizeof(double *) * data_depth);
    for (int i = 0; i < data_depth; ++i) network[i] = (double *) malloc(sizeof(double) * 10);
    for (int i = 0; i < data_depth; ++i) {
        for (int j = 0; j < 10; ++j) {
            double sum = 0;
            for (int k = 0; k < data_depth; ++k) {
                sum += result_cofactor_transpose[i][k] * result_matrix_ATB[k][j];
            }
            network[i][j] = sum;
        }
        printf("subprocess : process [ %d / %d ] \n", i + 1, data_depth);
    }
    printf("main process : ended multiple ( reverse ( A_transpose * A ), A * b ).\n");

    printf("main process : clean up memory used...\n");
    for (int i = 0; i < data_depth; ++i) free(result_cofactor_transpose[i]);
    free(result_cofactor_transpose);
    for (int i = 0; i < data_depth; ++i) free(result_matrix_ATB[i]);
    free(result_matrix_ATB);

    return network;

    /*
    // load or save cofactor transpose.
    double **result_cofactor_transpose = (double **) malloc(sizeof(double *) * data_depth);
    for (int i = 0; i < data_depth; ++i) result_cofactor_transpose[i] = (double *) malloc(sizeof(double) * data_depth);
    FILE *cofactor = fopen("cofactor_transpose", "rb");
    for (int i = 0; i < data_depth; ++i) fread(result_cofactor_transpose[i], 8, data_depth, cofactor);
    fclose(cofactor);


    // directly find pseudo inverse of A and save it.
    // multiple (ATA)-1AT
    double **result_pseudo_inverse = (double **) malloc(sizeof(double *) * data_depth);
    for (int i = 0; i < data_depth; ++i)
        result_pseudo_inverse[i] = (double *) malloc(sizeof(double) * train_data->data_count);
    for (int i = 0; i < data_depth; ++i) {
        for (int j = 0; j < train_data->data_count; ++j) {
            double sum = 0;
            for (int k = 0; k < data_depth; ++k) {
                sum += result_cofactor_transpose[i][k] * ((unsigned) (unsigned char) train_data->data[j][k]);
            }
            result_pseudo_inverse[i][j] = sum;
        }
        printf("subprocess : process [ %d / %d ] \n", i + 1, data_depth);
    }

    // load or save pseudo inverse.
    double **result_pseudo_inverse = (double **) malloc(sizeof(double *) * data_depth);
    for (int i = 0; i < data_depth; ++i)
        result_pseudo_inverse[i] = (double *) malloc(sizeof(double) * train_data->data_count);
    FILE *pseu = fopen("pseu", "rb");
    for (int i = 0; i < data_depth; ++i) fread(result_pseudo_inverse[i], 8, train_data->data_count, pseu);
    fclose(pseu);

    // process pseudo inverse times b.
    double **network = (double **) malloc(sizeof(double *) * data_depth);
    for (int i = 0; i < data_depth; ++i) network[i] = (double *) malloc(sizeof(double) * 10);
    for (int i = 0; i < data_depth; ++i) {
        for (int j = 0; j < 10; ++j) {
            double sum = 0;
            for (int k = 0; k < train_data->data_count; ++k) {
                sum += result_pseudo_inverse[i][k] * ((unsigned) (unsigned char) label_data->data[k][j]);
            }
            network[i][j] = sum;
        }
    }

    // free memory.
    for (int i = 0; i < data_depth; ++i) free(result_pseudo_inverse[i]);
    free(result_pseudo_inverse);

*/
}

void print_data_by_id(data_structure *train_data, label_structure *train_label, int select_id) {
    for (int j = 0; j < train_data->height; ++j) {
        for (int k = 0; k < train_data->length; ++k) {
            printf("%c",
                   ((unsigned) (unsigned char) train_data->data[select_id][j * train_data->length + k]) / 26 + '0');
        }
        printf("\n");
    }
    printf("data label = [");
    for (int j = 0; j < 10; ++j) {
        printf("%c", '0' + ((train_label->data[select_id][j]) ? 1 : 0));
    }
    printf("]\n\n");
}

int judge_matrix[10][10] = {0};

void test_validation(data_structure *test_data, label_structure *test_label, double **network) {
    memset(judge_matrix, 0, 4 * 100);
    unsigned data_depth = test_data->length * test_data->height;
    int revel_count = 0;
    for (int i = 0; i < test_data->data_count; ++i) {
        double result[10] = {0};
        for (int j = 0; j < 10; ++j) {
            double sum = 0;
            for (int k = 0; k < data_depth; ++k)
                sum += ((unsigned) (unsigned char) test_data->data[i][k]) * network[k][j];
            result[j] = sum;
        }
        int select = 0, target = 0;
        for (int j = 0; j < 10; ++j) {
            if (result[j] > result[select]) select = j;
            if (test_label->data[i][j]) target = j;
        }
        judge_matrix[target][select]++;
        if (select == target) revel_count++;
    }
    double acc = 100.00 * (double) revel_count / test_data->data_count;
    printf("[ %d / %d ] validation , accuracy = %lf %\n",
           revel_count, test_data->data_count, acc);
    printf("judge matrix : \n\t");
    for (int i = 0; i < 10; ++i) {
        printf("%d\t", i);
    }
    printf("\n");
    for (int i = 0; i < 10; ++i) {
        printf("%d\t", i);
        for (int j = 0; j < 10; ++j) {
            printf("%d\t", judge_matrix[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

// [30, 1, 32, 25, 25, 34, 30, 26, 28, 27]

int main() {
    // for trainning
    data_structure *train_data = load_data("train_file");
    label_structure *train_label = load_label("train_label");
    unsigned data_depth = train_data->height * train_data->length;

    double **net = calculate_network(train_data, train_label);

    FILE *network_record = fopen("network", "wb");
    for (int i = 0; i < data_depth; ++i) fwrite(net[i], 8, 10, network_record);
    fclose(network_record);

    for (int i = 0; i < data_depth; ++i) free(net[i]);
    free(net);


    // for testing
    data_structure *test_data = load_data("test_file");
    label_structure *test_label = load_label("test_label");
    // unsigned data_depth = test_data->height * test_data->length;

    double **neto = (double **) malloc(sizeof(double *) * data_depth);
    FILE *network_read = fopen("network", "rb");
    for (int i = 0; i < data_depth; ++i)
        neto[i] = (double *) malloc(sizeof(double) * 10), fread(neto[i], 8, 10, network_read);
    fclose(network_read);

    test_validation(test_data, test_label, neto);

    for (int i = 0; i < data_depth; ++i) free(neto[i]);
    free(neto);

    return 0;
}