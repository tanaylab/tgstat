#define BASE_ENV_FILE
#include "port.h"

void base_crash() {
#if DBG_ON
        if(stack_mode.to_act())
                stack_print(error_mode.out());
#endif // DBG_ON
        abort();
}
