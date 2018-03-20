#ifndef TGLUtil_port_h
#define TGLUtil_port_h 1

#ifndef TRACE
# define  TRACE    1
#endif
#ifndef DBG_ON
# define DBG_ON 0
#endif
#ifndef FUNCS
# define FUNCS     DBG_ON
#endif // FUNCS

#include <ctype.h>
#include <fcntl.h>
#include <float.h>
#include <iomanip>
#include <limits>
#include <limits.h>
#include <math.h>
#include <new>
#include <signal.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <unistd.h>
#include <sys/stat.h>
#include <time.h>

using namespace std;

#include "config.h"
#include "modes.h"
#include "error.h"
#include "init.h"
#include "trace.h"
#include "dbg.h"
#include "util.h"

#endif //  TGLUtil_port_h
