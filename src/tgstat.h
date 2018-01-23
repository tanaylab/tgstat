#ifndef TGSTAT_H_
#define TGSTAT_H_

#include <string>
#include <set>
#include <vector>
#include <pthread.h>
#include <semaphore.h>
#include <signal.h>
#include <time.h>

#include <R.h>
#include <Rinternals.h>

#include "Thread.h"

#ifdef length
#undef length
#endif
#ifdef error
#undef error
#endif

#include "TGLException.h"

using namespace std;

// should be used instead of R_CheckUserInterrupt. Throws exception if the command is interrupted.
void check_interrupt() throw (TGLException);

// adds timeout to the time that is already in req
void set_abs_timeout(int64_t delay_msec, struct timespec &req);

// sets timeout to req
void set_rel_timeout(int64_t delay_msec, struct timespec &req);

// returns true if current time exceeds start_time + delay
bool is_time_elapsed(int64_t delay_msec, const struct timespec &start_time);

// use rerror/verror instead of error!
void rerror(const char *fmt, ...);

void verror(const char *fmt, ...);

void vdebug(const char *fmt, ...);

// Use rprotect instead of PROTECT!
SEXP rprotect(SEXP &expr);

// Unprotect the last "count" object
void runprotect(int count);

// Unprotects object expr and sets it to R_NilValue. Works slower than runprotect(unsigned)!
void runprotect(SEXP &expr);

// Unprotects objects exprs and sets them to R_NilValue. Works slower than runprotect(unsigned)!
void runprotect(vector<SEXP> &exprs);

// Call runprotect_all if you wish to unprotect all object that are still protected
void runprotect_all();

inline bool is_R_var_char(char c) { return isalnum(c) || c == '_' || c == '.'; }

// the result is already protected
SEXP eval_in_R(SEXP parsed_command, SEXP envir);

// the result is already protected
SEXP run_in_R(const char *command, SEXP envir);

// This function writes R object to a file.
// Unlike R_Serialize function that just stops the execution if anything goes wrong (meaning: no clean up, destructors, etc.)
// RSaneSerialize throws an exception in case of error.
void RSaneSerialize(SEXP rexp, FILE *fp);
void RSaneSerialize(SEXP rexp, const char *fname);

// This function reads R object from a file. Object is expected to be saved using R's serialize() function or RSaneSerialize().
// The returned value is already protected.
// Unlike R_Unserialize function that just stops the execution if anything goes wrong (meaning: no clean up, destructors, etc.)
// RSaneUnserialize throws an exception in case of error.
SEXP RSaneUnserialize(FILE *fp);
SEXP RSaneUnserialize(const char *fname);

SEXP get_rvector_col(SEXP v, const char *colname, const char *varname, bool error_if_missing);

string get_bound_colname(const char *str, unsigned maxlen = 40);

template<typename T> void pack_data(void *&ptr, const T &data, size_t n) {
	size_t size = sizeof(data) * n;
	memcpy(ptr, &data, size);
	ptr = (char *)ptr + size;
}

template<typename T> void unpack_data(void *&ptr, T &data, size_t n) {
	size_t size = sizeof(data) * n;
	memcpy(&data, ptr, size);
	ptr = (char *)ptr + size;
}


#define MAX_KIDS 1000
#define rreturn(retv) { if (TGStat::is_kid()) exit(0); return(retv); }

// Define TGStat instance in your main function that is called by R.
// TGStat should be defined inside "try-catch" statement that catches TGLException.
// TGStat performs the following actions:
//   1. Installs a new SIGINT handler. ONE MUST CALL check_interrupt() INSTEAD OF R_CheckUserInterrupt()!!!!!!!
//   2. Installs out-of-memory handler.
//   3. Supresses the default error report behaviour.
//   4. Makes sure all file descriptors are closed on exit / error / interrupt.
//   5. Makes sure all objects are destructed on exit / error / interrupt.

class TGStat {
public:
	TGStat(SEXP _env);
	~TGStat();

	SEXP env() const { return m_env; }

    // Verifies that the data size does not exceed the maximum allowed.
    void verify_max_data_size(uint64_t data_size, const char *data_name = "Result");

    // true if debug prints are allowed
    bool debug() const { return m_debug; }

    // Returns the upper limit for data size
	uint64_t max_data_size() const { return m_max_data_size; }

    uint64_t rnd_seed() const { return m_rnd_seed; }
    void rnd_seed(uint64_t seed);

    static void set_alarm(int msecs);   // time is given in milliseconds
    static void reset_alarm();
    static int alarm_fired() { return s_sigalrm_fired; }

    static void prepare4multitasking();
    static pid_t launch_process();

    // returns false if one or more child processes have ended or true if the timeout has elapsed
    static bool wait_for_kid(int millisecs);

    // returns false if all the child processes have ended or true if the timeout has elapsed
    static bool wait_for_kids(int millisecs);

    // returns number of bytes read or 0 for EOF; the parent process that uses fifo does not need to call then wait_for_kids()
    static int read_multitask_fifo(void *buf, size_t bytes);
    static void write_multitask_fifo(const void *buf, size_t bytes);

    static bool is_kid() { return s_is_kid; }

    static void itr_idx(uint64_t idx) { s_shm->itr_idx[s_kid_index] = idx; }

    static uint64_t itr_idx_sum();   // sum of itr_idx over all the kids

    static int num_kids() { return s_kid_index; }

    static int num_kids_running() { return s_running_pids.size(); }

    static sem_t *shm_sem() { return s_shm_sem; }

protected:
    struct Shm {
        char          error_msg[10000];
        uint64_t      itr_idx[MAX_KIDS];          // used for progress report
    };

    struct SigBlocker {
        SigBlocker() {
            sigemptyset(&sigset);
            sigaddset(&sigset, SIGCHLD);
            sigaddset(&sigset, SIGINT);
            sigprocmask(SIG_BLOCK, &sigset, &oldsigset);
        }

        ~SigBlocker() { sigprocmask(SIG_UNBLOCK, &sigset, NULL); }

        sigset_t sigset;
        sigset_t oldsigset;
    };

	static struct sigaction     s_old_sigint_act;
    static struct sigaction     s_old_sigalrm_act;
    static struct sigaction     s_old_sigchld_act;
	static int                  s_ref_count;
	static int                  s_sigint_fired;
    static bool                 s_sigalrm_fired;
	static unsigned             s_protect_counter;

    static bool                 s_is_kid;
    static pid_t                s_parent_pid;
    static sem_t               *s_shm_sem;
    static sem_t               *s_fifo_sem;
    static int                  s_kid_index;
    static vector<pid_t>        s_running_pids;
    static Shm                 *s_shm;
    static int                  s_fifo_fd;

	SEXP                        m_env;
	mode_t                      m_old_umask;
	TGLException::Error_handler m_old_error_handler;
	new_handler                 m_old_new_handler;
	unsigned                    m_old_protect_count;
	set<int>                    m_old_open_fds;

    uint64_t                    m_rnd_seed;
    bool                        m_debug;
	uint64_t                    m_max_data_size;

	void load_options();

    static string  get_shm_sem_name();
    static string  get_fifo_sem_name();
    static string  get_fifo_name();
    static void    handle_error(const char *msg);
	static void    out_of_memory();
	static void    sigint_handler(int);
    static void    sigalrm_handler(int);
    static void    sigchld_handler(int);
	static void    get_open_fds(set<int> &fds);
    static void    check_kids_state(bool ignore_errors);

	friend void check_interrupt() throw (TGLException);
	friend SEXP rprotect(SEXP &expr);
	friend void runprotect(int count);
	friend void runprotect(SEXP &expr);
	friend void runprotect(vector<SEXP> &exprs);
	friend void runprotect_all();
	friend void rerror(const char *fmt, ...);
	friend void verror(const char *fmt, ...);
};

extern TGStat *g_tgstat;


// ------------------------------- IMPLEMENTATION --------------------------------

inline void set_abs_timeout(int64_t delay_msec, struct timespec &req)
{
	req.tv_nsec += delay_msec * 1000000L;
	req.tv_sec += req.tv_nsec / 1000000000L;
	req.tv_nsec %= 1000000000L;
}

inline void set_rel_timeout(int64_t delay_msec, struct timespec &req)
{
	req.tv_sec = delay_msec / 1000;
	req.tv_nsec = (delay_msec - req.tv_sec * 1000) * 1000000L;
}

inline bool is_time_elapsed(int64_t delay_msec, const struct timespec &start_time)
{
	struct timespec t1 = start_time;
	struct timespec t2;
	set_abs_timeout(delay_msec, t1);
	clock_gettime(CLOCK_REALTIME, &t2);
	return t2.tv_sec > t1.tv_sec || t2.tv_sec == t1.tv_sec && t2.tv_nsec > t1.tv_nsec;
}

inline uint64_t TGStat::itr_idx_sum()
{
    uint64_t res = 0;
    for (int i = 0; i < s_kid_index; ++i) 
        res += s_shm->itr_idx[i];
    return res;
}

#endif

