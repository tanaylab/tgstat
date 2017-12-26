#include "port.h"
BASE_CC_FILE
#include "BitVec.h"
#include "BitVecIter.h"

const uint ds_bitvec::tail_mask[] = {
                ~0U, 0x00000001U, 0x00000003U, 0x00000007U,
        0x0000000fU, 0x0000001fU, 0x0000003fU, 0x0000007fU,
        0x000000ffU, 0x000001ffU, 0x000003ffU, 0x000007ffU,
        0x00000fffU, 0x00001fffU, 0x00003fffU, 0x00007fffU,
        0x0000ffffU, 0x0001ffffU, 0x0003ffffU, 0x0007ffffU,
        0x000fffffU, 0x001fffffU, 0x003fffffU, 0x007fffffU,
        0x00ffffffU, 0x01ffffffU, 0x03ffffffU, 0x07ffffffU,
        0x0fffffffU, 0x1fffffffU, 0x3fffffffU, 0x7fffffffU
};
const uint1 ds_bitvec::ones[] = {
0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4, 
1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 
1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 
2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 
1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 
2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 
2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 
3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 
1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 
2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 
2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 
3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 
2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 
3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 
3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 
4, 5, 5, 6, 5, 6, 6, 7, 5, 6, 6, 7, 6, 7, 7, 8 
};


#define DBG_ASSERTEqualSizes \
        DBG_ASSERT(size_ == other.size_, \
                    "Binary op on sets with different universe sizes.")

ds_bitvec::ds_bitvec(int init_size) :
        size_(init_size),
        on_bits_(0),
        bits((init_size + IntBits - 1) / IntBits, uint(0))
{
//        ASSERT(size_ % IntBits <= int(numof(tail_mask)),
//                           "Private tail_mask array is not large enough");
}
ds_bitvec::ds_bitvec(const ds_bitvec &other, int offset, int init_size) :
        size_(init_size),
        on_bits_(-1),
        bits((init_size + IntBits - 1) / IntBits, uint(0))
{

        restrict(other, offset);
}
ds_bitvec::ds_bitvec(int init_size, int offset, const ds_bitvec &other) :
        size_(init_size),
        on_bits_(other.on_bits_),
        bits((init_size + IntBits - 1) / IntBits, uint(0))
{

        embed(offset, other);
}
/*ds_bitvec::ds_bitvec(istream &in) :
        size_(0),
        on_bits_(0),
        bits(0)
{

        in >> (*this);
}*/
void
ds_bitvec::set_ind(int ind) {

        int o = int_offset(ind);
        int i = int_index(ind);

        unsigned mask = (1 << o);

        if(!(bits[i] & mask)) {
                if(on_bits_ != -1)
                        on_bits_++;
                bits[i] |= mask;
        }
}
void
ds_bitvec::clr_ind(int ind) {

        int o = int_offset(ind);
        int i = int_index(ind);

        unsigned mask = (1 << o);

        if(bits[i] & mask) {
                if(on_bits_ != -1)
                        on_bits_--;
                bits[i] ^= mask;
        }
}
bool
ds_bitvec::flip_ind(int ind) {

        int o = int_offset(ind);
        int i = int_index(ind);

        unsigned mask = (1 << o);
        bool set_to = !(bits[i] & mask);
        bits[i] ^= mask;

        if(on_bits_ != -1) {
                if (set_to)
                        on_bits_++;
                else
                        on_bits_--;
        }

        return(set_to);
}
bool
ds_bitvec::get_ind(int ind) const {

        // (
        ASSERT(0 <= ind && ind < size_,
                "Bit index " << ind << " out of bounds [0.." << size_ << ")");
        // ]
        return((bits[int_index(ind)] >> int_offset(ind)) & 1);
}
void
ds_bitvec::set_size(int new_size) {

        if(new_size == size_)
                return;

//        DBG_ASSERT(new_size % IntBits <= int(numof(tail_mask)),
//                            "Private tail_mask array is not large enough");

        int old_bits_size = int(bits.size());
        int new_bits_size = (new_size + IntBits - 1) / IntBits;

        if(old_bits_size != new_bits_size) {
                bits.resize(new_bits_size);
	}

        int old_size = size_;
        size_ = new_size;
        if(new_size < old_size) {
                truncate_tail();
                on_bits_ = -1;
        }
}
ds_bitvec &ds_bitvec::set(bool bit) {

        if(bit) {
		fill_n(bits.begin(), bits.size(), ~0);
                truncate_tail();
                on_bits_ = size_;
        } else {
		fill_n(bits.begin(), bits.size(), 0);
                on_bits_ = 0;
        }

        return(*this);
}
ds_bitvec &ds_bitvec::clear() {

	fill_n(bits.begin(), bits.size(), 0);
        on_bits_ = 0;
        return(*this);
}
int ds_bitvec::on_bits() const {

        if(on_bits_ == -1) {
                int c = 0;
		for(vector<uint>::const_iterator i = bits.begin(); 
		    i != bits.end(); 
		    i++) {
			const uint1 *seg =(const uint1*)&(*i);
                        c += ones[seg[0]];
                        c += ones[seg[1]];
                        c += ones[seg[2]];
                        c += ones[seg[3]];
		}
                on_bits_ = c;
        }

        return(on_bits_);
}
int ds_bitvec::on_bits_and(const ds_bitvec &mask) const {

	int c = 0;
	vector<uint>::const_iterator i2 = mask.bits.begin(); 
	for(vector<uint>::const_iterator i1 = bits.begin(); 
	    i1 != bits.end(); 
	    i1++) {
		if(!(*i1 & *i2)) {
			i2++;
			continue;
		}
		const uint1 *seg1 =(const uint1*)&(*i1);
		const uint1 *seg2 =(const uint1*)&(*i2);
		c += ones[seg1[0]&seg2[0]];
		c += ones[seg1[1]&seg2[1]];
		c += ones[seg1[2]&seg2[2]];
		c += ones[seg1[3]&seg2[3]];
		i2++;
	}

        return(c);
}
int ds_bitvec::on_bits_and_pointers(list<pair<int, uint1> > &ps) const 
{
	int c = 0;
	uint1 *base = (uint1 *)&(*(bits.begin()));
	
	for(list<pair<int, uint1> >::iterator p = ps.begin();
	    p != ps.end();
	    p++) {
		c += ones[base[p->first]&p->second];
	}
        return(c);
}
int ds_bitvec::on_bits_minus(const ds_bitvec &mask) const {

	int c = 0;
	vector<uint>::const_iterator i2 = mask.bits.begin(); 
	for(vector<uint>::const_iterator i1 = bits.begin(); 
	    i1 != bits.end(); 
	    i1++) {
		if(!(*i1 & ~(*i2))) {
			i2++;
			continue;
		}
		const uint1 *seg1 =(const uint1*)&(*i1);
		const uint1 *seg2 =(const uint1*)&(*i2);
		c += ones[seg1[0]&~(seg2[0])];
		c += ones[seg1[1]&~(seg2[1])];
		c += ones[seg1[2]&~(seg2[2])];
		c += ones[seg1[3]&~(seg2[3])];
		i2++;
	}

        return(c);
}
ds_bitvec &ds_bitvec::noteq() {

        for(int i = 0; i < int(bits.size()); i++)
                bits[i] = ~(bits[i]);
        truncate_tail();

        if(on_bits_ != -1)
                on_bits_ = size_ - on_bits();

        return(*this);
}
ds_bitvec &ds_bitvec::restrict(const ds_bitvec &other, int offset) {

        DBG_ASSERT(offset + size_ <= other.size_,
                    "restricted range out of bounds:" << endl <<
                    "offset = " << offset <<
                    " size_ = " << size_ <<
                    " other.size_ = " << other.size_);

        int bit_shift = offset % IntBits;
        int j = offset / IntBits;

        if(bit_shift == 0) {
                for (int i = 0; i < int(bits.size()); i++, j++)
                        bits[i] = other.bits[j];
        } else {
                int rev_bit_shift = IntBits - bit_shift;
                int i = 0;
                for(; i < int(bits.size()) - 1; i++, j++) {
                        bits[i] = (other.bits[j] >> bit_shift);
                        bits[i] |= other.bits[j + 1] << rev_bit_shift;
                }

                bits[i] = (other.bits[j] >> bit_shift);
                if (j + 1 < int(other.bits.size()))
                        bits[i] |= other.bits[j + 1] << rev_bit_shift;
        }

        on_bits_ = -1;

        truncate_tail();

        return(*this);
}
ds_bitvec &ds_bitvec::embed(int offset, const ds_bitvec &other) {

        DBG_ASSERT(offset + other.size_ <= size_,
                    "restricted range out of bounds:" << endl <<
                    "offset = " << offset <<
                    " size_ = " << size_ <<
                    " other.size_ = " << other.size_);

	fill_n(bits.begin(), bits.size(), 0);

        int bit_shift = offset % IntBits;
        int j = offset / IntBits;

        if(bit_shift == 0) {
                for (int i = 0; i < int(other.bits.size()); i++, j++)
                        bits[j] = other.bits[i];
        } else {
                int rev_bit_shift = IntBits - bit_shift;
                int i = 0;
                for(; i < int(other.bits.size()) - 1; i++, j++) {
                        bits[j] |= (other.bits[i] << bit_shift);
                        bits[j + 1] = other.bits[i] >> rev_bit_shift;
                }

                bits[j] |= (other.bits[i] << bit_shift);
                if (j + 1 < int(bits.size()))
                        bits[j + 1] = other.bits[i] >> rev_bit_shift;
        }

        on_bits_ = other.on_bits_;

        truncate_tail();

        return(*this);
}
ds_bitvec &ds_bitvec::andeq(const ds_bitvec &other) {

        DBG_ASSERTEqualSizes;

        for(int i = 0; i < int(bits.size()); i++) {
                bits[i] &= other.bits[i];
	}
        on_bits_ = -1;
        return(*this);
}
ds_bitvec &ds_bitvec::oreq(const ds_bitvec &other) {

        DBG_ASSERTEqualSizes;

        for(int i = 0; i < int(bits.size()); i++)
                bits[i] |= other.bits[i];
        on_bits_ = -1;
        return(*this);
}
ds_bitvec &ds_bitvec::xoreq(const ds_bitvec &other) {

        DBG_ASSERTEqualSizes;

        for(int i = 0; i < int(bits.size()); i++)
                bits[i] ^= other.bits[i];

        on_bits_ = -1;
        return(*this);
}
ds_bitvec &ds_bitvec::minus(const ds_bitvec &other) {

        DBG_ASSERTEqualSizes;

        for(int i = 0; i < int(bits.size()); i++)
                bits[i] &= ~other.bits[i];

        on_bits_ = -1;
        return(*this);
}
ds_bitvec &
ds_bitvec::nand(const ds_bitvec &other) {

        DBG_ASSERTEqualSizes;

        for(int i = 0; i < int(bits.size()); i++) 
                bits[i] = ~(bits[i] & other.bits[i]);
        truncate_tail();

        on_bits_ = -1;
        return(*this);
}
ds_bitvec &
ds_bitvec::nor(const ds_bitvec &other) {

        DBG_ASSERTEqualSizes;

        for(int i = 0; i < int(bits.size()); i++) 
                bits[i] = ~(bits[i] | other.bits[i]);
        truncate_tail();

        on_bits_ = -1;
        return(*this);
}
ds_bitvec &
ds_bitvec::nxor(const ds_bitvec &other) {

        DBG_ASSERTEqualSizes;

        for(int i = 0; i < int(bits.size()); i++) 
                bits[i] = ~(bits[i] ^ other.bits[i]);
        truncate_tail();

        on_bits_ = -1;
        return(*this);
}
ds_bitvec &ds_bitvec::nminus(const ds_bitvec &other) 
{

        DBG_ASSERTEqualSizes;

        for(int i = 0; i < int(bits.size()); i++) 
                bits[i] = ~bits[i] | other.bits[i];
        truncate_tail();

        on_bits_ = -1;
        return(*this);
}
ds_bitvec &ds_bitvec::copy(const ds_bitvec &other) 
{
        
        DBG_ASSERTEqualSizes;

        size_ = other.size_;
        on_bits_ = other.on_bits_;
        bits = other.bits;

        return(*this);
}
ds_bitvec &ds_bitvec::ncopy(const ds_bitvec &other) 
{

        DBG_ASSERTEqualSizes;

        bits.resize(other.bits.size());
        for(int i = 0; i < int(bits.size()); i++)
                bits[i] = ~other.bits[i];
        truncate_tail();

        if(other.on_bits_ != -1)
                on_bits_ = size_ - other.on_bits_;
        else
                on_bits_ = -1;

        size_ = other.size_;
        return(*this);
}
ds_bitvec &ds_bitvec::sub(const ds_bitvec &other) 
{

        DBG_ASSERTEqualSizes;

        for(int i = 0; i < int(bits.size()); i++)
                bits[i] = ~bits[i] & other.bits[i];

        on_bits_ = -1;
        return(*this);
}
ds_bitvec &ds_bitvec::nsub(const ds_bitvec &other) 
{

        DBG_ASSERTEqualSizes;

        for(int i = 0; i < int(bits.size()); i++)
                bits[i] |= ~other.bits[i];
        truncate_tail();

        on_bits_ = -1;
        return(*this);
}
bool ds_bitvec::operator<=(const ds_bitvec &other) const
{

        DBG_ASSERTEqualSizes;

        for(int i = 0; i < int(bits.size()); i++)
                if((bits[i] & ~other.bits[i]) != 0)
                        return(false);
        return(true);
}
bool ds_bitvec::operator>=(const ds_bitvec &other) const
{

        return(other <= (*this));
}
bool ds_bitvec::operator==(const ds_bitvec &other) const
{
        DBG_ASSERTEqualSizes;

        for(int i = 0; i < int(bits.size()); i++)
                if(other.bits[i] != bits[i])
                        return(false);

        return(true);
}

size_t ds_bitvec::hash() const 
{
	size_t sum = 0;
	std::hash<unsigned int> h;
	for(vector<uint>::const_iterator i = bits.begin(); 
	    i != bits.end(); 
	    i++) {
		sum+= h(*i);
	}
	return(sum);
}
bool equal_to<const ds_bitvec *>::operator()(const ds_bitvec *const &b1, const ds_bitvec *const &b2) const {
	return((*b1) == (*b2));
}

size_t std::hash<const ds_bitvec *>::operator()(const ds_bitvec * const &b) const {
	return(b->hash());
}
size_t std::hash<ds_bitvec *>::operator()(ds_bitvec * const &b) const {
	return(b->hash());
}

bool equal_to<ds_bitvec>::operator()(ds_bitvec const &b1, ds_bitvec const &b2) const {
	return(b1 == b2);
}

size_t std::hash<ds_bitvec>::operator()(ds_bitvec const &b) const {
	return(b.hash());
}
ostream &operator<<(ostream &out, const ds_bitvec &bv) {
	if(bv.on_bits()) {
		for(ds_bitvec_iter i(bv); i.is_valid(); ++i) {
			out << int(i) << " ";
		}
	} else {
		out << "empty ";
	}
	return(out);
}
