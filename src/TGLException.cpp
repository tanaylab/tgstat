#include "TGLException.h"

// ------------------- TGLThrow ------------------

void TGLError(const char *format, ...)
{
	va_list ap;
	va_start(ap, format);
	TGLException e(-1, ap, format);
	va_end(ap);
	TGLException::s_error_handler(e);
}

void TGLError(int errcode, const char *format, ...)
{
	va_list ap;
	va_start(ap, format);
	TGLException e(errcode, ap, format);
	va_end(ap);
	TGLException::s_error_handler(e);
}


// ------------------- TGLException --------------

TGLException::Error_handler TGLException::s_error_handler = TGLException::throw_error_handler;

//void TGLException::base_error_handler(TGLException &e)
//{
//    BaseError(e.msg());
//}
//
void TGLException::throw_error_handler(TGLException &e)
{
	throw(e);
}

TGLException::Error_handler TGLException::set_error_handler(TGLException::Error_handler error_handler)
{
	Error_handler old_handler = s_error_handler;
	s_error_handler = error_handler;
	return old_handler;
}

TGLException::TGLException(int errcode, va_list &ap, const char *format) :
	m_errcode(errcode),
	m_type(typeid(TGLException::Unknown))
{
	msg(ap, format);
}

TGLException::TGLException(int errcode, const type_info &type, va_list &ap, const char *format) :
	m_errcode(errcode),
	m_type(type)
{
	msg(ap, format);
}

void TGLException::msg(va_list &ap, const char *format)
{
	char buf[MAX_ERROR_MSG_LEN + 1];
	vsnprintf(buf, MAX_ERROR_MSG_LEN + 1, format, ap);
	buf[MAX_ERROR_MSG_LEN] = 0;
	m_errorstr = buf;
}
