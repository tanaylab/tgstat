#include "port.h"
BASE_CC_FILE
#include "BitVecIter.h"

ds_bitvec_iter::ds_bitvec_iter(const ds_bitvec &vec) :
        bitvec(vec),
        valid(false)
{
        first(0);
}
ds_bitvec_iter::ds_bitvec_iter(const ds_bitvec &vec, int pos) :
        bitvec(vec),
        valid(false)
{
        first(pos);
}
int ds_bitvec_iter::first(int pos)
{
	if(pos >= bitvec.size()) {
		valid = 0;
		return(-1);
	}
	if(pos == 0) {
		current_int = 0;
		current_bit = 0;
		valid = !!bitvec.size_;

		if(valid && !(bitvec.bits[0] & 0x0001))
			++(*this);
	} else {
		seek(pos);
                if(!bitvec[pos]) {
			++(*this);
		}
	}

        return(*this);
}
int ds_bitvec_iter::last()
{

        valid = true;

        if(bitvec.size() % bitvec.IntBits == 0) {
                current_int = bitvec.bits.size();
                current_bit = 0;
        } else {
                current_int = bitvec.bits.size() - 1;
                current_bit = bitvec.size() - current_int * bitvec.IntBits;
        }
        return(--(*this));
}
ds_bitvec_iter::operator int() const {

        if(valid)
                return(current_int * bitvec.IntBits + current_bit);
        else
                return(-1);
}
void ds_bitvec_iter::seek(int to) {

        if(to == -1) {
                valid = false;
                return;
        }

        DBG_ASSERT(0 <= to && to < bitvec.size(),
                "Attempt to seek to position " << to <<
                " for a bit vector of size " << bitvec.size());

        current_int = to / bitvec.IntBits;
        current_bit = to % bitvec.IntBits;
        valid = true;
}
static int find_set_bit(uint bits) {

        int i = 0;
        while(!(bits & 0x0001)) {
                bits >>= 1;
                i++;
        }
        return i;
}

int ds_bitvec_iter::operator++() {

        ASSERT(valid, "Invalid iterator value at ++ operator");

        uint remaining_bits = 0;

        current_bit++;
        bool next_int = (current_bit == bitvec.IntBits);
        
        if(!next_int) {
                remaining_bits = bitvec.bits[current_int] >> current_bit;
                next_int = (remaining_bits == 0);
        }

        if(next_int) {
                do {
                        current_int++;
                        valid = (current_int < int(bitvec.bits.size()));
                } while(valid && bitvec.bits[current_int] == 0);

                if(valid)
                        current_bit = find_set_bit(bitvec.bits[current_int]);
        } else
                current_bit += find_set_bit(remaining_bits);

        if(valid) {
                int current = current_int * bitvec.IntBits + current_bit;
                return(current);
        } else
                return(-1);
}
int ds_bitvec_iter::find_set_bit_rev(uint bits)
{

        int i = bitvec.IntBits - 1;
        uint mask = (1 << i);

        while(!(bits & mask)) {
                bits <<= 1;
                i--;
        }

        return i;
}

int ds_bitvec_iter::operator--()
{

        ASSERT(valid, "Invalid iterator value at -- operator");

        uint remaining_bits = 0;

        bool prev_int = (current_bit == 0);
        if(!prev_int) {
                remaining_bits = bitvec.bits[current_int] <<
                                        (bitvec.IntBits - current_bit);
                prev_int = (remaining_bits == 0);
        }

        if(prev_int) {
                do {
                        current_int--;
                        valid = (current_int >= 0);
                } while (valid && bitvec.bits[current_int] == 0);
                if(valid)
                        current_bit =
                                find_set_bit_rev(bitvec.bits[current_int]);
        } else
                current_bit -= (bitvec.IntBits -
                                find_set_bit_rev(remaining_bits));

        if(valid) {
                int current = current_int * bitvec.IntBits + current_bit;
                return(current);
        } else
                return(-1);
}
