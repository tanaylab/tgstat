#ifndef __THREAD_H_INCLUDED__
#define __THREAD_H_INCLUDED__

#include <pthread.h>

struct MutexLocker {
	MutexLocker(pthread_mutex_t &mutex) : m_mutex(&mutex) { pthread_mutex_lock(m_mutex); }
	~MutexLocker() { pthread_mutex_unlock(m_mutex); }

	pthread_mutex_t *m_mutex;
};

struct SemLocker {
	SemLocker(sem_t *sem) : m_sem(sem) { sem_wait(m_sem); }
	~SemLocker() { sem_post(m_sem); }

	sem_t *m_sem;
};

#endif

