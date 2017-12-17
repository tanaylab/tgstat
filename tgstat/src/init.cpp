#include "port.h"

class EnvActions {
public:
        void (*funcs[64])();
        uint nfuncs;
public:
        ~EnvActions() {
                for(uint f = 0; f < nfuncs; f++)
                        (*funcs[f])();
        }
};
static EnvActions base_actions;
static ActiveMode base_error_mode("BASE_ERR_ACT", 1,
                                "BASE_ERR_LOG", 1,
                                "BASE_ERR_OUT", cerr);
#if DBG_ON
static ActiveMode base_stack_mode("BASE_STK_ACT", 1,
                                "BASE_STK_LOG", 0,
                                "BASE_STK_OUT", cerr);
#endif // DBG_ON

#if TRACE
static ActiveMode base_trace_mode("BASE_TRC_ACT", 1,
                                "BASE_TRC_LOG", 1,
                                "BASE_TRC_OUT", cerr);
#else // TRACE
InactiveMode base_trace_mode;
#endif // TRACE
#if DBG_ON
ActiveMode base_dbg_mode("BASE_PRN_ACT", 1,
                                "BASE_PRN_LOG", 1,
                                "BASE_PRN_OUT", cerr);
#else // DBG_ON
InactiveMode base_dbg_mode;
#endif // DBG_ON
ActiveMode &get_error_mode() { return(base_error_mode); }
#if DBG_ON
ActiveMode &get_stack_mode() { return(base_stack_mode); }
#endif // DBG_ON
#if TRACE
ActiveMode &get_trace_mode() { return(base_trace_mode); }
#else // TRACE
InactiveMode &get_trace_mode() { return(base_trace_mode); }
#endif // TRACE
#if DBG_ON
ActiveMode &get_dbg_mode() { return(base_dbg_mode); }
#else // DBG_ON
InactiveMode &get_dbg_mode() { return(base_dbg_mode); }
#endif // DBG_ON
void at_finish(void (*func)()) {
        ENV_FUNC("at_finish")
        /*if(base_actions.nfuncs >= numof(base_actions.funcs)) {
                cerr << "Asked to do too many things at finish" << endl;
                base_crash();
        }*/
        base_actions.funcs[base_actions.nfuncs++] = func;
}

bool arg_version(int argc, char *argv[]) {
        if(argc == 2 && !strcmp(argv[1], "version")) {
                cout << argv[0] << ' ' << version() << endl;
                return(true);
        } else {
                return(false);
	}
}
