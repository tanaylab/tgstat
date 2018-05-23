#ifndef TGLException_h_INCLUDED
#define TGLException_h_INCLUDED

#include <typeinfo>
#include <stdarg.h>
#include <stdio.h>
#include <string>

using namespace std;

// Use these four functions to handle fatal errors. They will generate a TGLException object and then
// call the installed error handler (TGLException::set_error_handler()).
// Use the template version of TGLThrow to keep track of the object that generated the exception
// (this object can be retrieved by TGLException::type())

void TGLAssert(bool cond, const char *format, ...);

void TGLAssert(bool cond, int errcode, const char *format, ...);

void TGLError(const char *format, ...);

void TGLError(int errcode, const char *format, ...);

template <typename T>
void TGLError(const char *format, ...);

template <typename T>
void TGLError(int errcode, const char *format, ...);


// ------------------- TGLException --------------

class TGLException {
public:
	enum { MAX_ERROR_MSG_LEN = 10000 };

	const char           *msg() const { return m_errorstr.c_str(); }
	int                   code() const { return m_errcode; }
	const std::type_info &type() const { return m_type; }

	typedef void (*Error_handler)(TGLException &);

//	static void base_error_handler(TGLException &e);
	static void throw_error_handler(TGLException &e);  // default error handler: throws TGLException

	// installs a new error handler, returns the old one
	static Error_handler set_error_handler(Error_handler error_handler);

private:
	class Unknown {};

	std::string           m_errorstr;
	int                   m_errcode;
	const std::type_info &m_type;
	static Error_handler  s_error_handler;

	TGLException(int errcode, va_list &ap, const char *format);
	TGLException(int errcode, const std::type_info &type, va_list &ap, const char *format);

	void msg(va_list &ap, const char *format);

    friend void TGLAssert(bool cond, const char *str, ...);
    friend void TGLAssert(bool cond, int errcode, const char *str, ...);
	friend void TGLError(const char *str, ...);
	friend void TGLError(int errcode, const char *str, ...);
	template <typename T> friend void TGLError(const char *str, ...);
	template <typename T> friend void TGLError(int errcode, const char *str, ...);
};


// ------------- implementation ---------------------

inline void TGLAssert(bool cond, const char *format, ...)
{
    if (!cond) {
    	va_list ap;
    	va_start(ap, format);
    	TGLException e(-1, ap, format);
    	va_end(ap);
    	TGLException::s_error_handler(e);
    }
}

inline void TGLAssert(bool cond, int errcode, const char *format, ...)
{
    if (!cond) {
    	va_list ap;
    	va_start(ap, format);
    	TGLException e(errcode, ap, format);
    	va_end(ap);
    	TGLException::s_error_handler(e);
    }
}

template <typename T> void TGLError(const char *format, ...)
{
	va_list ap;
	va_start(ap, format);
	TGLException e(-1, typeid(T), ap, format);
	va_end(ap);
	TGLException::s_error_handler(e);
}

template <typename T> void TGLError(int errcode, const char *format, ...)
{
	va_list ap;
	va_start(ap, format);
	TGLException e(errcode, typeid(T), ap, format);
	va_end(ap);
	TGLException::s_error_handler(e);
}

#endif
