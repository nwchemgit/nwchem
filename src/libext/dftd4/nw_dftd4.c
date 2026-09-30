/* DFTD4 interface for NWChem.
 * Developed with assistance from OpenAI Codex
 * Implementation and numerical validation
 * were reviewed and tested by MPH */
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "dftd4.h"

extern void nw_dftd4_get_params(const char *, double *, int *);

static void copy_error(dftd4_error error, char *message, int size)
{
    if (size <= 0) return;
    memset(message, ' ', (size_t)size);
    if (error) {
        char buffer[512] = {0};
        int n = (int)sizeof(buffer);
        dftd4_get_error(error, buffer, &n);
        n = (int)strnlen(buffer, sizeof(buffer));
        if (n > size) n = size;
        memcpy(message, buffer, (size_t)n);
    }
}

void nw_dftd4_eval_(const int64_t *natoms, const int64_t *numbers,
                    const double *positions, const double *charge,
                    const char *method, const int64_t *method_length,
                    const int64_t *mode, double *energy, double *derivative,
                    double *parameters, int64_t *version, int64_t *status,
                    char *message,
                    size_t ignored_method_length, size_t ignored_message_length)
{
    dftd4_error error = NULL;
    dftd4_structure mol = NULL;
    dftd4_model model = NULL;
    dftd4_param param = NULL;
    char name[256];
    int *numbers32 = NULL;
    int i, n, param_status;

    (void)ignored_method_length;
    (void)ignored_message_length;
    *status = 1;
    *energy = 0.0;
    *version = (int64_t)dftd4_get_version();
    memset(message, ' ', ignored_message_length);
    memset(parameters, 0, 6 * sizeof(*parameters));

    n = *method_length;
    if (n < 1 || n >= (int)sizeof(name)) {
        const char *text = "invalid DFTD4 method name length";
        memcpy(message, text, strlen(text));
        return;
    }
    memcpy(name, method, (size_t)n);
    name[n] = '\0';
    if (*natoms < 1 || *natoms > INT32_MAX) {
        const char *text = "invalid NWChem atom count for DFTD4";
        memcpy(message, text, strlen(text));
        return;
    }
    numbers32 = malloc((size_t)(*natoms) * sizeof(*numbers32));
    if (!numbers32) return;
    for (i = 0; i < *natoms; ++i) {
        if (numbers[i] < 1 || numbers[i] > INT32_MAX) {
            const char *text = "invalid atomic number passed to DFTD4";
            memcpy(message, text, strlen(text));
            goto done;
        }
        numbers32[i] = (int)numbers[i];
    }

    error = dftd4_new_error();
    if (!error) goto done;
    mol = dftd4_new_structure(error, (int)*natoms, numbers32, positions, charge,
                              NULL, NULL);
    if (dftd4_check_error(error) || !mol) goto done;
    model = dftd4_new_d4_model(error, mol);
    if (dftd4_check_error(error) || !model) goto done;
    param = dftd4_load_rational_damping(error, name, true);
    if (dftd4_check_error(error) || !param) goto done;
    param_status = 1;
    nw_dftd4_get_params(name, parameters, &param_status);
    if (param_status != 0) {
        const char *text = "could not read loaded DFTD4 parameters";
        memcpy(message, text, strlen(text));
        goto done;
    }

    if (*mode == 0) {
        dftd4_get_dispersion(error, mol, model, param, energy, NULL, NULL);
    } else if (*mode == 1) {
        dftd4_get_dispersion(error, mol, model, param, energy, derivative,
                             NULL);
    } else if (*mode == 2) {
        dftd4_get_numerical_hessian(error, mol, model, param, derivative);
    } else {
        const char *text = "invalid NWChem DFTD4 evaluation mode";
        memcpy(message, text, strlen(text));
        goto done;
    }
    if (!dftd4_check_error(error)) *status = 0;

done:
    if (*status != 0 && error) copy_error(error, message,
                                          ignored_message_length);
    dftd4_delete_param(&param);
    dftd4_delete_model(&model);
    dftd4_delete_structure(&mol);
    dftd4_delete_error(&error);
    free(numbers32);
}
