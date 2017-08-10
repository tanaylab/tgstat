#ifndef util_bitveciter_h
#define util_bitveciter_h 1

//We provide two almost identical iterators, one for the on bits and one for
//its complement. We are not sharing code to optimize performence (avoid
//virtuality and enable inlining). This iterator should provide one of
//the most possibly efficient bitvec implementations.

#include "BitVec.h"

class ds_bitvec_iter {
private:
        const ds_bitvec &bitvec;
        bool valid;
        int current_int;
        int current_bit;
public:
        ds_bitvec_iter(const ds_bitvec &vec);
        ds_bitvec_iter(const ds_bitvec &vec, int pos);
        int first(int pos = 0);
        int last();
        bool is_valid() const { return(valid); }
        operator int() const;
        void seek(int to);

        int operator++();
        int operator--();

private:
        int find_set_bit_rev(uint bits);
};
class ds_bitvec_notiter {

private:
        const ds_bitvec &bitvec;
        bool valid;
        int current_int;
        int current_bit;
public:
        ds_bitvec_notiter(const ds_bitvec &vec);
        int first();
        int last();
        bool is_valid() const { return(valid); }
        operator int() const;
        void seek(int to);

        int operator++();
        int operator--();

private:
        int find_unset_bit_rev(uint bits);
};

#endif // util_bitveciter_h
