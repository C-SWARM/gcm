#ifndef H__H__CONSTITUTIVE_MODEL_HANDLE__H__H
#define H__H__CONSTITUTIVE_MODEL_HANDLE__H__H

struct MATERIAL_ELASTICITY;
#ifndef TYPE_MATERIAL_ELASTICITY
#define TYPE_MATERIAL_ELASTICITY
typedef struct MATERIAL_ELASTICITY MATERIAL_ELASTICITY;
#endif

struct MATERIAL_CRYSTAL_PLASTICITY;
#ifndef TYPE_MATERIAL_CRYSTAL_PLASTICITY
#define TYPE_MATERIAL_CRYSTAL_PLASTICITY
typedef struct MATERIAL_CRYSTAL_PLASTICITY MATERIAL_CRYSTAL_PLASTICITY;
#endif

struct MATERIAL_CONTINUUM_DAMAGE;
#ifndef TYPE_MATERIAL_CONTINUUM_DAMAGEY
#define TYPE_MATERIAL_CONTINUUM_DAMAGE
typedef struct MATERIAL_CONTINUUM_DAMAGE MATERIAL_CONTINUUM_DAMAGE;
#endif

struct CRYSTAL_PLASTICITY_SOLVER_INFO;
#ifndef TYPE_CRYSTAL_PLASTICITY_SOLVER_INFO
#define TYPE_CRYSTAL_PLASTICITY_SOLVER_INFO
typedef struct CRYSTAL_PLASTICITY_SOLVER_INFO CRYSTAL_PLASTICITY_SOLVER_INFO;
#endif

struct CRYSTAL_PLASTICITY_SOLVER_INFO;
#ifndef TYPE_CRYSTAL_PLASTICITY_SOLVER_INFO
#define TYPE_CRYSTAL_PLASTICITY_SOLVER_INFO
typedef struct CRYSTAL_PLASTICITY_SOLVER_INFO CRYSTAL_PLASTICITY_SOLVER_INFO;
#endif

struct ELASTICITY;
#ifndef TYPE_ELASTICITY
#define TYPE_ELASTICITY
typedef struct ELASTICITY ELASTICITY;
#endif


struct CONTINUUM_DAMAGE;
#ifndef TYPE_CONTINUUM_DAMAGE
#define TYPE_CONTINUUM_DAMAGE
typedef struct CONTINUUM_DAMAGE CONTINUUM_DAMAGE;
#endif

typedef struct {
  MATERIAL_ELASTICITY            *mat_e;
  MATERIAL_CRYSTAL_PLASTICITY    *mat_p;
  MATERIAL_CONTINUUM_DAMAGE      *mat_d;
  CRYSTAL_PLASTICITY_SOLVER_INFO *solver_info;
  ELASTICITY                     *elasticity;
  CONTINUUM_DAMAGE               *damage;        
} CONSTITUTIVE_MODEL_PACKS;
#endif
