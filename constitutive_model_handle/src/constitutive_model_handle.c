#include <stdio.h>
#include "constitutive_model_handle.h"

int constitutive_model_handle_init(CONSTITUTIVE_MODEL_PACKS *cm_pack)
{
  int err = 0;
  cm_pack->mat_e       = NULL;
  cm_pack->mat_p       = NULL;
  cm_pack->mat_d       = NULL;
  cm_pack->solver_info = NULL;
  cm_pack->elasticity  = NULL;
  cm_pack->damage      = NULL;
  return err;
}
