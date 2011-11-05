#include <stdio.h>

#define LINE 512

/* Read in graphical data from a file, list the "marker" lines",
   and echo the desired data on standard output */
main()
{
	FILE *strm;
	char filen[80], buf[LINE];
	char marks[100][10], wantmark[100][10];
	int markno, nplot;
	int i,j;
	double xval, yval;

	fprintf(stderr,"Enter filename\n");
	scanf("%s", filen);
	strm = fopen(filen, "r");
	markno = 0;
	while (fgets(buf, LINE, strm) != NULL) {
		if (sscanf(buf, " marker %s", marks[markno]) == 1) {
			fprintf(stderr, "Marker: %s\n", marks[markno]);
			markno++;
		}
	}
	fclose(strm);
	strm = fopen(filen, "r"); /* Primitive rewinding */
	fprintf(stderr, "How many of these do you want plotted?\n");
	scanf("%d", &nplot);
	fprintf(stderr, "Enter them, in order of appearance in file\n");
	j = 0;
	for (i=0; i < nplot; i++) {
	    scanf("%s", wantmark[i]);
	    printf("#\n");
	    for (; j < markno; j++) {
		if (!strcmp(marks[j],wantmark[i])) {
		    while (fgets(buf, LINE, strm) != NULL) {
			if (sscanf(buf, "%lf %lf", &xval, &yval) == 2) {
				printf("%g %g\n", xval, yval);
			} else {
				j++;
				break;
			}
		    }
		    break;
		} else {
		    while (fgets(buf, LINE, strm) != NULL) {
			if (sscanf(buf, "%lf %lf", &xval, &yval) != 2)
				break;
		    }
		}
	    }
	}
}
/* $Id$ */
