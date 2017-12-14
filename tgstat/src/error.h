#ifndef base_error_h
#define base_error_h 1

// Basic assertion
#define BaseError(Message) \
        if(error_mode.to_log()) { \
		error_mode.out() << HERE << Message << endl; \
                base_crash(); \
        } else if(always()) { \
                base_crash(); \
        } else  // Eats the trailing ';'

#define ASSERT(Cond, Message) \
        if (!(Cond)) {\
		BaseError(Message); \
        } else  // Eats the trailing ';'


#if DBG_ON
extern ActiveMode &get_dbg_mode();
#define DBG_ASSERT(Cond, Message) \
        if(dbg_mode.to_act() && !(Cond)) \
                BaseError(Message); \
        else    // Eats the trailing ;
#else // DBG_ON
extern InactiveMode &get_dbg_mode();
#define DBG_ASSERT(Cond, Message)
#endif // DBG_ON

#define BasePrn BaseLog(dbg_mode)

extern ActiveMode &get_error_mode();

extern void base_crash();

#endif // base_error_h
