#ifndef base_config_h
#define base_config_h 1

typedef signed char int1;
typedef unsigned char uint1;
typedef signed short int2;
typedef unsigned short uint2;
typedef signed long int4;
typedef unsigned long uint4;
typedef unsigned int uint;
typedef float real;
typedef float real4;
typedef double real8;
#define _INT(ATTR)          INT_##ATTR
#define _UINT(ATTR)         UINT_##ATTR
#define _INT1(ATTR)         SCHAR_##ATTR
#define _UINT1(ATTR)        UCHAR_##ATTR
#define _INT2(ATTR)         SHRT_##ATTR
#define _UINT2(ATTR)        USHRT_##ATTR
#define ENENVNT4(ATTR)         LONG_##ATTR
#define _UINT4(ATTR)        ULONG_##ATTR
#define _REAL(ATTR)         FLT_##ATTR
#define _REAL4(ATTR)        FLT_##ATTR
#define _REAL8(ATTR)        DBL_##ATTR
#define FLT_SAME               FLT_EPSILON
#define DBL_SAME               FLT_EPSILON
#define always() true
#define never() false
#if !ENV_HAS_FOR_SCOPE
# define for    if(never()) ; else for
#endif // ENV_HAS_FOR_SCOPE
#define ENV_NO_ARGS
#if ENV_HAS_INST_TEMPLATE

#define ENV_INST_START(FUNC_ID, FUNC_ARGS)
#define ENV_INST_CLASS(CLASS, ID, ARGS)        template class CLASS;
#define ENV_INST_FINISH

#else // ENV_HAS_INST_TEMPLATE

#define ENV_INST_START(FUNC_ID, FUNC_ARGS) \
                        void FUNC_ID FUNC_ARGS { SEA_FUNC(#FUNC_ID)
#define ENV_INST_CLASS(CLASS, ID, ARGS) CLASS ID ARGS; BaseUseVar(ID);
#define ENV_INST_FINISH }

#endif // ENV_HAS_INST_TEMPLATE
#define BaseUseVar(V)   (void)((void*)(&V))

#endif // base_config_h
