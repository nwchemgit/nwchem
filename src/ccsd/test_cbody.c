#include <stdio.h>
#include <stdlib.h>

int main(int argc, char* argv[])
{
    int ncor, nocc, nvir;

    if (argc!=3) {
        printf("Usage: ./test_cbody nocc nvir\n");
        return argc;
    } else {
        ncor = 0;
        nocc = atoi(argv[1]);
        nvir = atoi(argv[2]);
    }

    if (nocc<1 || nvir<1) {
        printf("Arguments must be non-negative!\n");
        return 1;
    }

    printf("Test driver for cbody with nocc=%d, nvir=%d\n", nocc, nvir);

    return 0;
}
