#ifndef base_modes_h
#define base_modes_h 1

class ActiveMode {
public:
        bool to_act() const;
        void act();
        void deact();
        bool to_log() const;
        void log();
        void nolog();
        ostream &out() const;
        ActiveMode(const char *act_env, int act_def,
                        const char *log_env, int log_def,
                        const char *out_env, ostream &out_def);
        ~ActiveMode();
private:
        ActiveMode(const ActiveMode &);
        void operator=(const ActiveMode &);
        int act_;
        int log_;
        ostream *out_;
};
class InactiveMode {
public:
        bool to_act() const;
        void act();
        void deact();
        bool to_log() const;
        void log();
        void nolog();
        ostream &out() const;
};
#define BaseLog(Mode) if(!(Mode).to_log()) { } else (Mode).out()

inline bool ActiveMode::to_act() const { return(act_ > 0); }
inline void ActiveMode::act() { act_++; }
inline void ActiveMode::deact() { act_--; }
inline bool ActiveMode::to_log() const { return(log_ > 0); }
inline void ActiveMode::log() { log_++; }
inline void ActiveMode::nolog() { log_--; }
inline ostream &ActiveMode::out() const { return(*out_); }
inline bool InactiveMode::to_act() const { return(false); }
inline void InactiveMode::act() { }
inline void InactiveMode::deact() { }
inline bool InactiveMode::to_log() const { return(false); }
inline void InactiveMode::log() { }
inline void InactiveMode::nolog() { }
inline ostream &InactiveMode::out() const { return(cerr); }

#endif // base_modes_h
