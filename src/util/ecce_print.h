/*
 $Id$
 */

void ecce_print_module_entry(const char *module);
void ecce_print_module_exit(const char *module, const char *status);
void ecce_print1(const char *keyword, int ma_type,
		 const void *data, int dim1);
void ecce_print2(const char *keyword, int ma_type, 
		 const void *data, int ld1, int dim1, int dim2);
void ecce_print_control(int new, int *old);
void ecce_print_file_open(const char *filename);
void ecce_print_file_close(void);
void ecce_print_echo_input(const char *filename);
void ecce_print_echo_string(const char *string);
void ecce_print_version(const char *string);
logical is_ecce_print_on(void);
