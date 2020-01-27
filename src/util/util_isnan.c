/* $Id$ */
/* function to detect NaNs by calling C isnan */
#include <stdio.h>
#include <math.h>
#include "typesf2c.h"

Logical FATR
util_isnan_(DoublePrecision *number_in) {
  int yesno = (isnan(*number_in)) ? 1: 0;
  return yesno;
}
