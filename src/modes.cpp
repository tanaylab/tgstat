#include "port.h"

ActiveMode::ActiveMode(
        const char *act_env, int act_def,
        const char *log_env, int log_def,
        const char *out_env, ostream &out_def
){
        const char *cp;

        cp = getenv(act_env);
        act_ = (cp ? atoi(cp) : act_def);

        cp = getenv(log_env);
        log_ = (cp ? atoi(cp) : log_def);

        cp = getenv(out_env);
        if(!cp) {
                out_ = &out_def;
        } else if(!strcmp(cp, "out")
               || !strcmp(cp, "OUT")) {
                out_ = &cout;
        } else if(!strcmp(cp, "err")
               || !strcmp(cp, "ERR")) {
                out_ = &cerr;
        } else {
                out_ = new ofstream(cp);
                if(!out_) {
                        out_ = &cerr;
                        cerr << "ERROR: Can't open file " << cp
                                << " for " << out_env << endl;
                        exit(1);
                }
        }
}

ActiveMode::~ActiveMode() {
        if(out_ != &cout && out_ != &cerr)
                delete(out_);
}
ActiveMode::ActiveMode(const ActiveMode &) {
        cerr << "PANIC: ActiveMode::ActiveMode called" << endl;
        exit(EXIT_FAILURE);
}

void ActiveMode::operator=(const ActiveMode &) {
        cerr << "PANIC: ActiveMode::operator= called" << endl;
        exit(EXIT_FAILURE);
}
