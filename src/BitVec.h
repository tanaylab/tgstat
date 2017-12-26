#ifndef util_bitvec_h
#define util_bitvec_h 1

#include <limits.h>
#include <vector>
#include <list>
#include <unordered_map>

class ds_bitvec {

        friend class ds_bitvec_iter;
        friend class ds_bitvec_noiter;
private:
        #ifdef WORD_BIT
        enum { IntBits = WORD_BIT };
        #else // WORD_BIT
        enum { IntBits = CHAR_BIT * sizeof(int) };
        #endif // WORD_BIT
        class ds_bit {
        private:
                int ind;
                ds_bitvec *This;
        public:
                ds_bit(ds_bitvec *init_this, int init_ind) :
			ind(init_ind),
			This(init_this)
		{
		}
                operator bool() const {
			return(This->get_ind(ind));
		}
                bool operator=(bool bit) {
			if(bit)
				This->set_ind(ind);
			else
				This->clr_ind(ind);
			return(bit);
		}
                bool operator=(const ds_bitvec::ds_bit &bit){
			return ((*this) = bool(bit));
		}
                bool operator++(int) {
			This->set_ind(ind);
			return(true);
		}
                bool operator--(int) {
			This->clr_ind(ind);
			return(false);
		}
                bool flip() {
			return(This->flip_ind(ind));
		}
        };
        friend class ds_bit;
public:
        int size() const { return(size_); }
	uint get_word(int i) {
		return(bits[i]);
	}
	uint1 get_char(int i) {
		uint1 *base = (uint1 *)&(*(bits.begin()));
		return(base[i]);
	}
        ds_bit operator[](int ind) {
		return(ds_bit(this, ind)); 
	}
        bool operator[](int ind) const {
		return(get_ind(ind)); 
	}
        ds_bitvec(int init_size = 0);
        ds_bitvec(const ds_bitvec &other, int offset, int init_size);
        ds_bitvec(int init_size, int offset, const ds_bitvec &other);
        //ds_bitvec(istream &in);
        void set_size(int new_size);

        ds_bitvec &set(bool bit = true);
        ds_bitvec &clear();
        int on_bits() const;
        int on_bits_and(const ds_bitvec &mask) const;
	int on_bits_and_pointers(list<pair<int, uint1> > &ps) const;
        int on_bits_minus(const ds_bitvec &mask) const;
        ds_bitvec &noteq();
        operator const void*() const {
		return(on_bits() ? this : 0); 
	}
        ds_bitvec &restrict(const ds_bitvec &other, int offset);

        ds_bitvec &embed(int offset, const ds_bitvec &other);

        ds_bitvec &andeq(const ds_bitvec &other);
        ds_bitvec &oreq(const ds_bitvec &other);
        ds_bitvec &xoreq(const ds_bitvec &other);
        ds_bitvec &minus(const ds_bitvec &other);
        ds_bitvec &nand(const ds_bitvec &other);
        ds_bitvec &nor(const ds_bitvec &other);
        ds_bitvec &nxor(const ds_bitvec &other);
        ds_bitvec &nminus(const ds_bitvec &other);
        ds_bitvec &copy(const ds_bitvec &other);
        ds_bitvec &ncopy(const ds_bitvec &other);
        ds_bitvec &sub(const ds_bitvec &other);
        ds_bitvec &nsub(const ds_bitvec &other);
        bool operator<=(const ds_bitvec &other) const;
        bool operator>=(const ds_bitvec &other) const;

        bool operator==(const ds_bitvec &other) const;
        bool operator!=(const ds_bitvec &other) const {
		return(!(*this == other));
	}

        size_t hash() const;
private:
        int size_;
        mutable int on_bits_;
	vector<uint> bits;

        static const uint tail_mask[];
        static const uint1 ones[];

        void set_ind(int ind);
        void clr_ind(int ind);
        bool flip_ind(int ind);
        bool get_ind(int ind) const;
        void truncate_tail() {
		if(size_)
			bits[bits.size() - 1] &= tail_mask[size_ % IntBits];
	}

        static int int_index(int i) {
		return(i / ds_bitvec::IntBits);
	}
        static int int_offset(int i) {
		return(i % ds_bitvec::IntBits);
	}
};

ostream &operator<<(ostream &out, const ds_bitvec &bv);

namespace std {
template<> class equal_to<ds_bitvec> {
public:
	bool operator()(ds_bitvec const &b1, ds_bitvec const &b2) const;
};
template<> class equal_to<const ds_bitvec *> {
public:
	bool operator()(const ds_bitvec * const &b1, 
				const ds_bitvec * const &b2) const;
};
}

//bool equal_to<const ds_bitvec *>::operator()(const ds_bitvec * const &b1, const ds_bitvec * const &b2) const;
namespace std {
template<> class hash<const ds_bitvec *> {
public:
	size_t operator()(const ds_bitvec * const &b1) const;
};
template<> class hash<ds_bitvec *> {
public:
	size_t operator()(ds_bitvec * const &b1) const;
};
template<> class hash<ds_bitvec> {
public:
	size_t operator()(ds_bitvec const &b1) const;
};
}

//size_t hash<const ds_bitvec *>::operator()(const ds_bitvec *const &b) const;

#endif // util_bitvec_h
