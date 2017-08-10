#ifndef base_stack_h
#define base_stack_h 1

#if DBG_ON 
extern ActiveMode & get_stack_mode();
extern void  stack_push(const char *file, int line, const char *func);
extern void  stack_pop();
extern void  stack_print(ostream &out);
class StackTrace {
public:
         StackTrace(const char *file, int line, const char *func);
         ~StackTrace();
};
#if FUNCS
# define ENV_FUNC(F)       \
        static const char *const FUNC_NAME = F; \
        UseVar(FUNC_NAME);
#else
# define ENV_FUNC(F)
# define FUNC_NAME ""
#endif
#define STK_TRACE(F) \
        static const char *const FUNC_NAME = F; \
        StackTrace stack_trace(__FILE__, __LINE__, FUNC_NAME);
#else // DBG_ON
#define ENV_FUNC(F)
#define FUNC_NAME ""
#define STK_TRACE(F)
#endif // DBG_ON

#if DBG_ON
inline StackTrace::StackTrace(const char *file, int line, const char *func) {
        if(stack_mode.to_act())
                stack_push(file, line, func);
}

inline StackTrace::~StackTrace() {
        if(stack_mode.to_act())
                stack_pop();
}
#endif // DBG_ON

extern const char *const  sea_msg();
#define ScmStk ScmLog(stack_mode)
#ifdef FUNC_NAME
# define HERE __FILE__ << '(' << __LINE__ << "): "
#else // FUNC_NAME
# define HERE __FILE__ << '(' << __LINE__ << "): " << FUNC_NAME << ": "
#endif // FUNC_NAME

#endif // base_stack_h
