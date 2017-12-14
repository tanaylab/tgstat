#ifndef base_util_h
#define base_util_h 1

inline double max(double f1, double f2) { return(f1 > f2 ? f1 : f2); }
inline float max(float f1, float f2) { return(f1 > f2 ? f1 : f2); }
inline int max(int f1, int f2) { return(f1 > f2 ? f1 : f2); }
inline float max(int f1, float f2) { return(f1 > f2 ? f1 : f2); }
inline float max(float f1, int f2) { return(f1 > f2 ? f1 : f2); }

inline double min(double f1, double f2) { return(f1 < f2 ? f1 : f2); }
inline float min(float f1, float f2) { return(f1 < f2 ? f1 : f2); }
inline int min(int f1, int f2) { return(f1 < f2 ? f1 : f2); }
inline float min(int f1, float f2) { return(f1 < f2 ? f1 : f2); }
inline float min(float f1, int f2) { return(f1 < f2 ? f1 : f2); }

inline void log_sum_log(float &l1, float l2) {
	if(l1 > l2) {
		if(!isinf(l2)) {
			l1 += log(1+ exp(l2-l1));
		}
	} else {
		if(isinf(l1)) {
			l1 = l2;
		} else {
			l1 = l2 + log(1 + exp(l1-l2));
		}
	}
}

inline void log_minus_log(float &l1, float l2) {
	if(1 < exp(l2 - l1)) {
		if(l2 - l1 < fabs(l1*0.0001)) {
			l1 = -_REAL(MAX);
				
		} else {
//			cerr << "NAN at log minus log " 
//			<< l1 << " " << l2 << " " << l2 - l1 << endl;
			l1 = -_REAL(MAX);
		}
	} else {
		l1 += log(1 - exp(l2-l1));
	}
}

inline float log_one_minus_exp(float l1) {
	return(log(1-exp(l1)));
}

inline float log_one_minus(float l1) {
	if(l1 < 0.001) {
		return(-l1);
	} else {
		return(log(1 - l1));
	}
}

inline void log_sum_log(double &l1, double l2) {
	if(l1 > l2) {
		l1 += log(1+ exp(l2-l1));
	} else {
		l1 = l2 + log(1 + exp(l1-l2));
	}
}

inline void log_minus_log(double &l1, double l2) {
	if(1 < exp(l2 - l1)) {
		cerr << "NAN at log minus log " 
			<< l1 << " " << l2 << " " << l2 - l1 << endl;
	}
	l1 += log(1 - exp(l2-l1));
}

inline double log_one_minus_exp(double l1) {
	return(log(1-exp(l1)));
}

inline double log_one_minus(double l1) {
	if(l1 < 0.001) {
		return(-l1);
	} else {
		return(log(1 - l1));
	}
}

template<class T, class T1, class T2>
class triplet {

public:
	T first;
	T1 second;
	T2 third;

	triplet(T o, T1 o1, T2 o2) :
		first(o),
		second(o1),
		third(o2)
	{}
};

#endif // base_uil_h

