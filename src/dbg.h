#ifndef base_dbg_h
#define base_dbg_h 1

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

#endif // base_dbg_h
