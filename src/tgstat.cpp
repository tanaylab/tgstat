#include <ctype.h>
#include <dirent.h>
#include <errno.h>
#include <fcntl.h>
#include <limits>
#include <new>
#include <setjmp.h>
#include <stdarg.h>
#include <unistd.h>
#include <cstdint>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/wait.h>

#ifndef _POSIX_C_SOURCE
	#define _POSIX_C_SOURCE 199309
	#include <time.h>
	#undef _POSIX_C_SOURCE
#endif

#if defined(__APPLE__)
    #include <libproc.h>
    #include <sys/proc_info.h>
#endif

#ifndef R_NO_REMAP
#  define R_NO_REMAP
#endif
#include <R.h>
#include <Rinternals.h>
#include <Rinterface.h>
#include <Rembedded.h>
#include <R_ext/Parse.h>

#ifdef length
#undef length
#endif

#ifdef error
#undef error
#endif

#include "tgstat.h"

using namespace std;

struct sigaction         TGStat::s_old_sigint_act;
struct sigaction         TGStat::s_old_sigalrm_act;
struct sigaction         TGStat::s_old_sigchld_act;
int                      TGStat::s_ref_count = 0;
int                      TGStat::s_sigint_fired = 0;
bool                     TGStat::s_sigalrm_fired = 0;
unsigned                 TGStat::s_protect_counter = 0;
bool                     TGStat::s_is_kid = false;
pid_t                    TGStat::s_parent_pid = 0;
sem_t                   *TGStat::s_shm_sem = SEM_FAILED;
sem_t                   *TGStat::s_fifo_sem = SEM_FAILED;
int                      TGStat::s_kid_index;
vector<pid_t>            TGStat::s_running_pids;
TGStat::Shm             *TGStat::s_shm = (TGStat::Shm *)MAP_FAILED;
int                      TGStat::s_fifo_fd = -1;

TGStat *g_tgstat = NULL;

TGStat::TGStat(SEXP _env) :
	m_env(_env)
{
// disable R check stack limit: required if eval is called not from the main thread
//R_CStackLimit=-1;
	if (!s_ref_count) {
        GetRNGstate();

		m_old_umask = umask(07);

        s_sigint_fired = 0;
        s_sigalrm_fired = 0;

        s_is_kid = false;
        s_kid_index = 0;
        s_parent_pid = getpid();
        s_shm_sem = SEM_FAILED;
        s_fifo_sem = SEM_FAILED;
        s_shm = (Shm *)MAP_FAILED;
        s_fifo_fd = -1;
        s_running_pids.clear();

		m_old_error_handler = TGLException::set_error_handler(TGLException::throw_error_handler);

		struct sigaction new_act;

		// install a new SIGINT handler
		new_act.sa_handler = sigint_handler;
		sigemptyset(&new_act.sa_mask);
		new_act.sa_flags = SA_RESTART;
		sigaction(SIGINT, &new_act, &s_old_sigint_act);

        // install a new SIGALRM handler
        new_act.sa_handler = sigalrm_handler;
        sigemptyset(&new_act.sa_mask);
        new_act.sa_flags = SA_RESTART;
        sigaction(SIGALRM, &new_act, &s_old_sigalrm_act);

        // install a new SIGCHLD handler
        new_act.sa_handler = sigchld_handler;
        sigemptyset(&new_act.sa_mask);
        new_act.sa_flags = SA_RESTART | SA_NOCLDSTOP;
        sigaction(SIGCHLD, &new_act, &s_old_sigchld_act);

		// record the currently opened file descriptors
		get_open_fds(m_old_open_fds);

		load_options();
	}

	s_ref_count++;

	// deal with PROTECT / UNPROTECT
	m_old_protect_count = s_protect_counter;

	if (s_ref_count == 1)
		g_tgstat = this;
}

TGStat::~TGStat()
{
	s_ref_count--;

	if (!s_ref_count) {
        // if this is a child, do not detach from shared memory and do not deallocate the semaphore:
        // if exception is thrown ~RdbInitializer is called first and then the child might need
        // to write the error into the shared memory
        if (!s_is_kid) {
            PutRNGstate();

            if (s_shm_sem != SEM_FAILED) {
                SemLocker sl(s_shm_sem);
                SigBlocker sb;

                // kill all the remaining child processes
                for (vector<pid_t>::const_iterator ipid = s_running_pids.begin(); ipid != s_running_pids.end(); ++ipid) {
                    vdebug("Forcefully terminating process %d\n", *ipid);
                    kill(*ipid, SIGTERM);
                }
            }

            // after SIGTERM is sent to all the kids let's wait till sigchld_hander() burries them all
            while (1) {
                SigBlocker sb;
                check_kids_state(true);
                if (s_running_pids.empty())
                    break;

                vdebug("Waiting for %ld child processes to end\n", s_running_pids.size());
                sigsuspend(&sb.oldsigset);
            }

            if (s_shm_sem != SEM_FAILED)
                sem_close(s_shm_sem); // semaphore should be already unlinked, only need to close it

            if (s_fifo_sem != SEM_FAILED)
                sem_close(s_fifo_sem);

            if (s_shm != (Shm *)MAP_FAILED)
                munmap((char *)s_shm, sizeof(Shm));   // needs to be char * for some versions of Solaris

            unlink(get_fifo_name().c_str());
        }

        if (s_fifo_fd != -1)
            close(s_fifo_fd);

        TGLException::set_error_handler(m_old_error_handler);

        // reset alarm, otherwise it might fire later when we exit from the library
        alarm(0);

		// install old signal handlers
		sigaction(SIGINT, &s_old_sigint_act, NULL);
        sigaction(SIGALRM, &s_old_sigalrm_act, NULL);
        sigaction(SIGCHLD, &s_old_sigchld_act, NULL);

		// close all file descriptors opened during the session
		set<int> open_fds;
		get_open_fds(open_fds);
		for (set<int>::const_iterator ifd = open_fds.begin(); ifd != open_fds.end(); ++ifd) {
			if (m_old_open_fds.find(*ifd) == m_old_open_fds.end())
				close(*ifd);
		}

		umask(m_old_umask);

		// do not revert to R's default error report
	}

	// deal with PROTECT / UNPROTECT
	runprotect(s_protect_counter - m_old_protect_count);
	s_protect_counter = m_old_protect_count;

	if (!s_ref_count)
		g_tgstat = NULL;
}

string TGStat::get_shm_sem_name()
{
	char buf[100];
	snprintf(buf, sizeof(buf), "/tgstat_shm_sem_%d", (int)getpid());
	return buf;
}

string TGStat::get_fifo_sem_name()
{
	char buf[100];
	snprintf(buf, sizeof(buf), "/tgstat_fifo_sem_%d", (int)getpid());
	return buf;
}

string TGStat::get_fifo_name()
{
	char buf[100];
    snprintf(buf, sizeof(buf), "/tmp/tgstat_fifo_%d", s_is_kid ? (int)getppid() : (int)getpid());
	return buf;
}

void TGStat::prepare4multitasking()
{
    vdebug("Cleaning old semaphores\n");
	if (s_shm_sem == SEM_FAILED) {
		sem_unlink(get_shm_sem_name().c_str()); // remove a semaphore if it was somehow not cleaned from the previous invocation of the lib
		if ((s_shm_sem = sem_open(get_shm_sem_name().c_str(), O_CREAT | O_EXCL, 0644, 1)) == SEM_FAILED)
			verror("sem_open failed: %s", strerror(errno));

		// Open a semaphore and right after that unlink it. The semaphore will be useful up until
		// the last process holding it (i.e. keeping it open) dies.
		sem_unlink(get_shm_sem_name().c_str());
	}

    if (s_fifo_sem == SEM_FAILED) {
        sem_unlink(get_fifo_sem_name().c_str()); // remove a semaphore if it was somehow not cleaned from the previous invocation of the lib
        if ((s_fifo_sem = sem_open(get_fifo_sem_name().c_str(), O_CREAT | O_EXCL, 0644, 1)) == SEM_FAILED)
            verror("sem_open failed: %s", strerror(errno));

        // Open a semaphore and right after that unlink it. The semaphore will be useful up until
        // the last process holding it (i.e. keeping it open) dies.
        sem_unlink(get_fifo_sem_name().c_str());
    }

    vdebug("Creating FIFO channel\n");
    if (s_fifo_fd == -1) {
        unlink(get_fifo_name().c_str());

        if (mkfifo(get_fifo_name().c_str(), 0666) == -1)
            verror("mkfifo of file %s failed: %s", get_fifo_name().c_str(), strerror(errno));

        if ((s_fifo_fd = open(get_fifo_name().c_str(), O_RDONLY | O_NONBLOCK)) == -1)
            verror("open of fifo %s for read failed: %s", get_fifo_name().c_str(), strerror(errno));

#ifdef F_SETPIPE_SZ
        fcntl(s_fifo_fd, F_SETPIPE_SZ, 1048576);   // set pipe size to 1 Mb
#endif
    }

    vdebug("Allocating shared memory for internal communication\n");
	if (s_shm == (Shm *)MAP_FAILED) {
		s_shm = (Shm *)mmap(NULL, sizeof(Shm), PROT_READ | PROT_WRITE, MAP_SHARED | MAP_ANONYMOUS, -1, 0);

		if (s_shm == (Shm *)MAP_FAILED)
            verror("Failed to allocate shared memory: %s", strerror(errno));

		s_shm->error_msg[0] = '\0';
		for (int i = 0; i < MAX_KIDS; ++i)
			s_shm->itr_idx[i] = 0;
	}
}

pid_t TGStat::launch_process()
{
	if (s_shm_sem == SEM_FAILED || s_fifo_sem == SEM_FAILED || s_shm == (Shm *)MAP_FAILED || s_fifo_fd == -1)
		verror("Not ready for multitasking");

	if (s_kid_index >= MAX_KIDS) 
		verror("Too many child processes");

    vdebug("SemLock\n");

	check_interrupt();

	{
		SemLocker sl(s_shm_sem);
		if (s_shm->error_msg[0])
			verror("%s", s_shm->error_msg);
	}

    vdebug("fork\n");
	pid_t pid = fork(); 

	if (pid == -1)
        verror("fork failed: %s", strerror(errno));

	if (pid) { // a parent process
        vdebug("%d: child process %d has been launched\n", getpid(), pid);
        s_running_pids.push_back(pid);
        ++s_kid_index;
    } else {   // a child process
		s_is_kid = true;

		sigaction(SIGINT, &s_old_sigint_act, NULL);
        sigaction(SIGALRM, &s_old_sigalrm_act, NULL);
		sigaction(SIGCHLD, &s_old_sigchld_act, NULL);

		SEXP r_multitasking_stdout = Rf_GetOption1(Rf_install("tgs_multitasking_stdout"));
        int devnull;

        if ((devnull = open("/dev/null", O_RDWR)) == -1)
            verror("Failed to open /dev/null");

        if (!Rf_isLogical(r_multitasking_stdout) || !(int)LOGICAL(r_multitasking_stdout)[0])
            dup2(devnull, STDOUT_FILENO);

        dup2(devnull, STDIN_FILENO);
        dup2(devnull, STDERR_FILENO);
        close(devnull);

        // fifo was open for read by the parent
        close(s_fifo_fd);

        if ((s_fifo_fd = open(get_fifo_name().c_str(), O_WRONLY)) == -1)
            verror("open of fifo %s for write failed: %s", get_fifo_name().c_str(), strerror(errno));
	}

	return pid;
}

void TGStat::check_kids_state(bool ignore_errors)
{
    int status;
    pid_t pid;

    while ((pid = waitpid((pid_t)-1, &status, WNOHANG)) > 0) {
        vdebug("pid %d has ended\n", pid);
        for (vector<pid_t>::iterator ipid = s_running_pids.begin(); ipid != s_running_pids.end(); ++ipid) {
            if (*ipid == pid) {
                vdebug("pid %d was identified as a child process\n", pid);
                swap(*ipid, s_running_pids.back());
                s_running_pids.pop_back();
                if (!ignore_errors && !WIFEXITED(status) && WIFSIGNALED(status) && WTERMSIG(status) != TGS_EXIT_SIG)
                    verror("Child process %d ended unexpectedly", (int)pid);
                break;
            }
        }
    }
}

bool TGStat::wait_for_kid(int millisecs)
{
    struct timespec timeout, remaining;
    set_rel_timeout(millisecs, timeout);

    while (1) {
        vdebug("SIGINT fired? %d\n", s_sigint_fired);
        check_interrupt();

        uint64_t num_running_pids = s_running_pids.size();
        check_kids_state(false);

        {
            SemLocker sl(s_shm_sem);
            if (s_shm->error_msg[0])
                verror("%s", s_shm->error_msg);
        }

        if (s_running_pids.empty() || num_running_pids > s_running_pids.size()) {
            vdebug("still running %ld child processes\n", s_running_pids.size());
            return false;
        }

        vdebug("still running %ld child processes (%d, ...)\n", s_running_pids.size(), s_running_pids.front());

        if (nanosleep(&timeout, &remaining))
            timeout = remaining;
        else
            break;
    }
    return true;
}

bool TGStat::wait_for_kids(int millisecs)
{
    struct timespec timeout, remaining;
    set_rel_timeout(millisecs, timeout);

    while (1) {
        vdebug("SIGINT fired? %d\n", s_sigint_fired);
        check_interrupt();
        check_kids_state(false);

        {
            SemLocker sl(s_shm_sem);
            if (s_shm->error_msg[0])
                verror("%s", s_shm->error_msg);
        }

        if (s_running_pids.empty()) {
            vdebug("No more running child processes\n");
            return false;
        }

        vdebug("still running %ld child processes (%d, ...)\n", s_running_pids.size(), s_running_pids.front());

        if (nanosleep(&timeout, &remaining))
            timeout = remaining;
        else
            break;
    }
    return true;
}

int TGStat::read_multitask_fifo(void *buf, uint64_t bytes)
{
    bool eof_reached = false;
    int retv;
    uint64_t readlen = 0;
    fd_set rfds;
    struct timeval tv;

    while (bytes > readlen) {
        tv.tv_sec = 1;  // wake up every second and check the interrupt
        tv.tv_usec = 0;

        FD_ZERO(&rfds);
        FD_SET(s_fifo_fd, &rfds);

        // remember: "select" is not restartable, i.e. it may be interrupted by a signal such as SIGCHLD even if we use SA_RESTART in the signal handler
        retv = select(s_fifo_fd + 1, &rfds, NULL, NULL, &tv);

        if (retv == -1) {
            if (errno != EINTR)
                verror("select on fifo failed: %s", strerror(errno));
        } else if (retv == 1) {
            retv = read(s_fifo_fd, buf, bytes - readlen);

            if (retv == -1) {
                if (errno != EAGAIN && errno != EWOULDBLOCK)
                    verror("read from fifo failed: %s", strerror(errno));
            } else {
                buf = (char *)buf + retv;
                readlen += retv;

                if (!retv)
                    eof_reached = true;
            }
        }

        check_interrupt();

        if (s_shm->error_msg[0]) {
            SemLocker sl(s_shm_sem);
            verror("%s", s_shm->error_msg);
        }

        check_kids_state(false);

        // If we exit on EOF without checking s_running_pids, we might miss the error message that the kid generated before it exited.
        // This happens because on error the kid first closes the FIFO and only then generates an error message and exits. 
        if (eof_reached && s_running_pids.empty())
            return readlen;
    }
    return readlen;
}

void TGStat::write_multitask_fifo(const void *buf, uint64_t bytes)
{
    SemLocker sl(s_fifo_sem);
    if (write(s_fifo_fd, buf, bytes) == -1)
        verror("write to fifo failed: %s", strerror(errno));
}

void TGStat::handle_error(const char *msg)
{
	if (s_is_kid) {
		{
			SemLocker sl(s_shm_sem);
			if (!s_shm->error_msg[0]) { // write an error message only if there were no error messages before
				strncpy(s_shm->error_msg, msg, sizeof(s_shm->error_msg) - 1);
				s_shm->error_msg[sizeof(s_shm->error_msg) - 1] = '\0';
			}
		}
		rexit();
	} else {
        Rf_errorcall(R_NilValue, "%s", msg);
    }
}

void TGStat::set_alarm(int msecs)
{
    struct itimerval timer;

    timer.it_interval.tv_sec = 0;
    timer.it_interval.tv_usec = 0;

    timer.it_value.tv_sec = msecs / 1000;
    timer.it_value.tv_usec = (msecs % 1000) * 1000;

    setitimer(ITIMER_REAL, &timer, NULL);
}

void TGStat::reset_alarm()
{
    s_sigalrm_fired = 0;

    struct itimerval timer;

    timer.it_interval.tv_sec = 0;
    timer.it_interval.tv_usec = 0;

    timer.it_value.tv_sec = 0;
    timer.it_value.tv_usec = 0;

    setitimer(ITIMER_REAL, &timer, NULL);
}

void TGStat::load_options()
{
	SEXP rvar;

    rvar = Rf_GetOption1(Rf_install("tgs_debug"));
    if (Rf_isLogical(rvar))
        m_debug = (int)LOGICAL(rvar)[0];
    else
        m_debug = false;

    int num_cores = max(1, (int)sysconf(_SC_NPROCESSORS_ONLN));
    rvar = Rf_GetOption1(Rf_install("tgs_max.processes"));
    if (Rf_xlength(rvar) && (Rf_isNumeric(rvar) || Rf_isInteger(rvar))) {
        m_num_processes = Rf_asInteger(rvar);
        m_num_processes = max(m_num_processes, 1);
        m_num_processes = min(m_num_processes, num_cores);
    } else
        m_num_processes = num_cores;
}

void TGStat::rnd_seed(uint64_t seed)
{
    SEXP e;
    SEXP sseed = Rf_install("set.seed");
    PROTECT(e = Rf_lang2(sseed, Rf_ScalarInteger(seed)));
    R_tryEval(e, m_env, NULL);
    UNPROTECT(1);
    GetRNGstate();
}

void TGStat::sigint_handler(int)
{
	++s_sigint_fired;

    // Normally this condition should be always true since the kid installs the default handler for SIGINT.
    // However due to race condition the old handler might still be in use.
    if (getpid() == s_parent_pid)
        REprintf("CTL-C!\n");
}

void TGStat::sigalrm_handler(int)
{
	s_sigalrm_fired = true;
}

void TGStat::sigchld_handler(int)
{
}

void TGStat::get_open_fds(set<int> &fds)
{
#if defined(__APPLE__)
    // This absolutely irrational code with all those funny reallocations and multiplication by 32 (haeh?) was inherited from here:
    //      https://opensource.apple.com/source/Libc/Libc-825.26/darwin/proc_listpidspath.c.auto.html
    // 
    // It would be much more logical to call proc_pidinfo twice: the first time to get buf size and the second time to get the list
    // of file descriptors. And indeed some internet sources advice to do that. However what is rational does not work and gives some
    // 240 open file descriptors most of them are not vnodes at all. And even some of those who are vnodes, are complete junk. Smells
    // like a memory leak. Probably the multiplication by 32 is really needed.
    //
    // All these problems come from the simple reason there is not manual for proc_pidinfo. Go figure why...
	int	buf_used;
    int fds_size = 0;
    int num_fds;
    unique_ptr<char[]> buf;

	// get list of open file descriptors
	buf_used = proc_pidinfo(getpid(), PROC_PIDLISTFDS, 0, NULL, 0);
	if (buf_used <= 0)
        return;

	while (1) {
		if (buf_used > fds_size) {
			// if we need to allocate [more] space
			while (buf_used > fds_size)
				fds_size += (sizeof(struct proc_fdinfo) * 32);

            buf = unique_ptr<char[]>(new char[fds_size]);
		}

		buf_used = proc_pidinfo(getpid(), PROC_PIDLISTFDS, 0, buf.get(), fds_size);
		if (buf_used <= 0)
            return;

		if ((buf_used + sizeof(struct proc_fdinfo)) >= fds_size) {
			// if not enough room in the buffer for an extra fd
			buf_used = fds_size + sizeof(struct proc_fdinfo);
			continue;
		}

		num_fds = buf_used / sizeof(struct proc_fdinfo);
		break;
	}

    struct proc_fdinfo *fdinfo = (struct proc_fdinfo *)buf.get();
    for (int i = 0; i < num_fds; ++i) {
        if (fdinfo[i].proc_fdtype == PROX_FDTYPE_VNODE)
            fds.insert(fdinfo[i].proc_fd);
    }
#else

#ifdef __sun
    #ifdef __XOPEN_OR_POSIX
        #define _dirfd(dir) (dir->d_fd)
    #else
        #define _dirfd(dir) (dir->dd_fd)
    #endif
#else
    #define _dirfd(dir) dirfd(dir)
#endif
	DIR *dir = opendir("/proc/self/fd");
	struct dirent *dirp;

	fds.clear();
    if (dir) {
    	while ((dirp = readdir(dir))) {
    		char *endptr;
    		int fd = strtol(dirp->d_name, &endptr, 10);
    		if (!*endptr && fd != _dirfd(dir)) // name is a number (it can be also ".", "..", whatever...)
    			fds.insert(fd);
    	}

        closedir(dir);
    }
#endif
}

void check_interrupt()
{
	if (TGStat::s_sigint_fired)
		TGLError("Command interrupted!");
}

void rerror(const char *fmt, ...)
{
	va_list ap;
	char buf[1000];

	va_start(ap, fmt);
	vsnprintf(buf, sizeof(buf), fmt, ap);
	va_end(ap);

	TGStat::handle_error(buf);
}

void verror(const char *fmt, ...)
{
	va_list ap;
	char buf[1000];

	va_start(ap, fmt);
	vsnprintf(buf, sizeof(buf), fmt, ap);
	va_end(ap);

	if (TGStat::s_ref_count)
		TGLError("%s", buf);
	else
		TGStat::handle_error(buf);
}

void vdebug(const char *fmt, ...)
{
    if (g_tgstat->debug()) {
        struct timeval tmnow;
        struct tm *tm;
        char buf[1000];

        gettimeofday(&tmnow, NULL);
        tm = localtime(&tmnow.tv_sec);
        strftime(buf, sizeof(buf), "%H:%M:%S", tm);
        if (TGStat::is_kid())
            REprintf("[DEBUG pid %d %s.%03d] ", getpid(), buf, (int)(tmnow.tv_usec / 1000));
        else
            REprintf("[DEBUG %s.%03d] ", buf, (int)(tmnow.tv_usec / 1000));

        va_list ap;
    	va_start(ap, fmt);
    	vsnprintf(buf, sizeof(buf), fmt, ap);
        va_end(ap);
        REprintf("%s", buf);

        if (!*fmt || (*fmt && fmt[strlen(fmt) - 1] != '\n'))
            REprintf("\n");
    }
}

SEXP rprotect(SEXP &expr)
{
	if (expr != R_NilValue) {
		TGStat::s_protect_counter++;
		return PROTECT(expr);
	}
	return expr;
}

void runprotect(unsigned count)
{
	if (TGStat::s_protect_counter < count)
		Rf_errorcall(R_NilValue, "Number of calls to unprotect exceeds the number of calls to protect\n");
	UNPROTECT(count);
	TGStat::s_protect_counter -= count;
}

void runprotect(SEXP &expr)
{
	if (expr != R_NilValue) {
		if (TGStat::s_protect_counter < 1)
			Rf_errorcall(R_NilValue, "Number of calls to unprotect exceeds the number of calls to protect\n");
		UNPROTECT_PTR(expr);
		expr = R_NilValue;
		TGStat::s_protect_counter--;
	}
}

void runprotect(vector<SEXP> &exprs)
{
	for (vector<SEXP>::iterator iexpr = exprs.begin(); iexpr != exprs.end(); ++iexpr)
		runprotect(*iexpr);
}

void runprotect_all()
{
	if (TGStat::s_protect_counter)
		UNPROTECT(TGStat::s_protect_counter);
	TGStat::s_protect_counter -= 0;
}

const char *get_groot(SEXP envir)
{
	// no need to protect the returned value
	SEXP groot = Rf_findVar(Rf_install("GROOT"), envir);

	if (!Rf_isString(groot))
		verror("GROOT variable does not exist");

	return CHAR(STRING_ELT(groot, 0));
}

const char *get_glib_dir(SEXP envir)
{
	// no need to protect the returned value
	SEXP glibdir = Rf_findVar(Rf_install(".GLIBDIR"), envir);

	if (!Rf_isString(glibdir))
		verror(".GLIBDIR variable does not exist");

	return CHAR(STRING_ELT(glibdir, 0));
}

SEXP eval_in_R(SEXP parsed_command, SEXP envir)
{
    int check_error;
    SEXP res;
    SEXP err_call;

    rprotect(res = R_tryEval(parsed_command, envir, &check_error));
    if (check_error)
    {
        // Get the last error message
        PROTECT(err_call = Rf_lang1(Rf_install("geterrmessage")));
        SEXP err_msg = R_tryEval(err_call, R_GlobalEnv, &check_error);
        UNPROTECT(1);

        if (!check_error && TYPEOF(err_msg) == STRSXP && LENGTH(err_msg) > 0)
        {
            verror(CHAR(STRING_ELT(err_msg, 0)));
        }
        else
        {
            verror("R evaluation error: Unknown error");
        }
    }
    return res;
}

SEXP run_in_R(const char *command, SEXP envir)
{
	SEXP expr;
	SEXP parsed_expr = R_NilValue;
    SEXPCleaner parsed_expr_cleaner(parsed_expr);
	ParseStatus status;

	rprotect(expr = RSaneAllocVector(STRSXP, 1));
	SET_STRING_ELT(expr, 0, Rf_mkChar(command));
	rprotect(parsed_expr = R_ParseVector(expr, -1, &status, R_NilValue));
	if (status != PARSE_OK)
		verror("Failed to parse expression \"%s\"", command);

	return eval_in_R(VECTOR_ELT(parsed_expr, 0), envir);
}

struct RSaneSerializeData {
	SEXP  rexp;
	FILE *fp;
};

static void RSaneSerializeCallback(void *_data)
{
	RSaneSerializeData *data = (RSaneSerializeData *)_data;
	struct R_outpstream_st out;
	R_InitFileOutPStream(&out, data->fp, R_pstream_xdr_format, 2, NULL, NULL);
	R_Serialize(data->rexp, &out);
}

void RSaneSerialize(SEXP rexp, FILE *fp)
{
	RSaneSerializeData data;

	data.rexp = rexp;
	data.fp = fp;
    Rboolean ok = R_ToplevelExec(RSaneSerializeCallback, &data);
	if (ok == FALSE)
		// We would like to print now the error contained in R_curErrorBuf(), however this error is automatically printed by R_ToplevelExec
		// and there's no way to prevent it without heavy hacking. On the other hand we want to abort the execution on error.
		// Solution: write a different error message. :)
		verror("Execution aborted");
}

void RSaneSerialize(SEXP rexp, const char *fname)
{
	FILE *fp = fopen(fname, "w");

	if (!fp)
		verror("Failed to open file %s: %s", fname, strerror(errno));

	RSaneSerialize(rexp, fp);
	fclose(fp);
}

struct RSaneUnserializeData {
	FILE *fp;
	SEXP  retv;
};

static void RSaneUnserializeCallback(void *_data)
{
	RSaneUnserializeData *data = (RSaneUnserializeData *)_data;
	struct R_inpstream_st in;
	R_InitFileInPStream(&in, data->fp, R_pstream_xdr_format, NULL, NULL);
	rprotect(data->retv = R_Unserialize(&in));
}

SEXP RSaneUnserialize(FILE *fp)
{
    RSaneUnserializeData data;

    data.fp = fp;
    data.retv = R_NilValue;

    Rboolean ok = R_ToplevelExec(RSaneUnserializeCallback, &data);
	if (ok == FALSE)
		// We would like to print now the error contained in R_curErrorBuf(), however this error is automatically printed by R_ToplevelExec
		// and there's no way to prevent it without heavy hacking. On the other hand we want to abort the execution on error.
		// Solution: write a different error message. :)
		verror("Execution aborted");
	runprotect(1);
	return data.retv;
}

SEXP RSaneUnserialize(const char *fname)
{
	FILE *fp = fopen(fname, "r");

	if (!fp)
		verror("Failed to open file %s: %s", fname, strerror(errno));

	SEXP retv = RSaneUnserialize(fp);

	fclose(fp);
	return retv;
}

struct RSaneAllocVectorData {
    SEXPTYPE type;
    R_xlen_t len;
    SEXP     retv;
};

static void RSaneAllocVectorCallback(void *_data)
{
	RSaneAllocVectorData *data = (RSaneAllocVectorData *)_data;
    data->retv = Rf_allocVector(data->type, data->len);
}

SEXP RSaneAllocVector(SEXPTYPE type, R_xlen_t len)
{
    RSaneAllocVectorData data;

    data.type = type;
    data.len = len;
    Rboolean ok = R_ToplevelExec(RSaneAllocVectorCallback, &data);
    if (!ok)
        verror("Allocation failed");
    return data.retv;
}

SEXP get_rvector_col(SEXP v, const char *colname, const char *varname, bool error_if_missing)
{
	SEXP colnames = Rf_getAttrib(v, R_NamesSymbol);

	if (!Rf_isVector(v) || (Rf_length(v) && (!Rf_isString(colnames) || Rf_length(colnames) != Rf_length(v))) || (!Rf_length(v) && !Rf_isNull(colnames)))
		verror("Invalid format of %s", varname);

	int numcols = Rf_isNull(colnames) ? 0 : Rf_length(colnames);

	for (int i = 0; i < numcols; i++) {
		if (!strcmp(CHAR(STRING_ELT(colnames, i)), colname))
			return VECTOR_ELT(v, i);
	}

	if (error_if_missing)
		verror("Invalid format of %s: missing %s column", varname, colname);
	return R_NilValue;
}

string get_bound_colname(const char *str, unsigned maxlen)
{
	string colname;

	maxlen = max(maxlen, 4u);
	if (strlen(str) > maxlen) {
		colname.assign(str, maxlen - 3);
		colname += "...";
	} else
		colname = str;
	return colname;
}

