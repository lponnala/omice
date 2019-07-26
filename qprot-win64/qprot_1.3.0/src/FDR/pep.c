#include "pep.h"
int nrow(FILE *fp) {
    char buf[10000];
    int n = 0;
    while(fgets(buf, sizeof(buf), fp) != NULL) n++;
    return n;
}

int newlinechar(char *buf, int k) {
    int i;
    int found = 0;
    for(i=0;i<k;i++) {
        if(buf[i] == '\n') {
            found = 1;
            break;
        }
    }   
    return found;
}

int ncol(FILE *fp) {
    char buf[10000];
    int i,cont = 0;
    fgets(buf, sizeof(buf), fp);
    for(i=0;i<10000;i++) {
        if(buf[i] == '\t') cont++;
        if(buf[i] == '\0') break;
    }
    return cont;
}

int main(int argc, char **argv) {

    /* Data */
    FILE *fpm;
    /* DATA ratio; */
    char outfile[100];

    if (argc != 2) {
        fprintf(stderr, "usage: getfdr <matrixData> \n");
        return 1;
    }

    FILE *fpi = fopen(argv[1], "r");
    np = nrow(fpi);
    rewind(fpi);
    nq = ncol(fpi)+1;
    fclose(fpi);

    fpm = fopen(argv[1], "r");

    np -= 1;
    PRINTDENSITY = 1;
    THEORETICAL = 1;  

    /* Read Data */
    read_matrix_data(fpm); 
    init_param();  
    fit_FDR();
    computeFDR();

    strcpy(outfile, argv[1]);
    strcat(outfile, "_fdr");
    FILE *fp_out = fopen(outfile, "w");
    print_result(fp_out);
    fclose(fpm);
    fclose(fp_out);

    if(PRINTDENSITY) {
        strcpy(outfile, argv[1]);
        strcat(outfile, "_density");
        FILE *fp_dens = fopen(outfile, "w");
        print_density(fp_dens);
        fclose(fp_dens);
    }
    return 0;

}





