#ifndef stdalg_ds_RaList_h
#define stdalg_ds_RaList_h 1

#include <vector>

template<class T> class RaList;
template<class T> class RaListIter;

template<class T> 
class RaListNode {

public:
	friend class RaList<T>;
	friend class RaListIter<T>;
private:

	T obj;
	int id;
	int next_id;
	int prev_id;
public:
	T &get_obj() {
		return(obj);
	}
	const T &get_obj() const {
		return(obj);
	}

	int get_id() {
		return(id);
	}
	RaListNode() :
		id(-1)
	{}
};

template<class T> class RaListIter;

template<class T>
class RaList {

	friend class RaListIter<T>;

public:
	typedef RaListIter<T> iterator;
//	typedef typename RaList<T>::iterator;

protected:

	vector<RaListNode<T> > m_nodes;

	int m_front_id;
	int m_back_id;

	uint m_size;

public:

	iterator begin();

	iterator end();

	uint get_size() const { return(m_size); }

	uint get_id_space_size() const { return(m_nodes.size()); }

	const T &operator[](int elem_id) const {
		return(m_nodes[elem_id].get_obj());
	}
	T &operator[](int elem_id) {
		return(m_nodes[elem_id].get_obj());
	}

	bool is_member(int elem_id) const {
		return(uint(elem_id) < m_nodes.size() 
			&& m_nodes[elem_id].id != -1);
	}

	int push_front(const T &obj, int elem_id = -1);
	int push_back(const T &obj, int elem_id = -1);
	int insert_after(const T &obj, int after_id, int elem_id);
	void remove(int elem_id);

	int get_front_id() const {
		return(m_front_id);
	}
	int get_back_id() const {
		return(m_back_id);
	}
	RaListNode<T> &get_front() {
		return(m_nodes[m_front_id]);
	}

	RaListNode<T> &get_back() {
		return(m_nodes[m_back_id]);
	}

	int get_prev_id(int elem_id) const {
		return(m_nodes[elem_id].prev_id);
	}
	int get_next_id(int elem_id) const {
		return(m_nodes[elem_id].next_id);
	}

	void reserve(int id_space_size) {
		m_nodes.resize(id_space_size);
	}

	void clear() {
		m_size = 0;
		m_front_id = m_back_id = -1;
		//Note the most efficient thing to do, it can be improved
		//by looking at size/id_space_size and removing one by one
		//in some cases. 
		fill(m_nodes.begin(), m_nodes.end(), RaListNode<T>());
	}

	//allocating id if needed
	RaList(int id_space_size = 0) :
		m_nodes(id_space_size),
		m_front_id(-1),
		m_back_id(-1),
		m_size(0)
	{
		cerr << "RaList with " << id_space_size << endl;
	}

private:
	int get_free_id();
};

#endif // stdalg_ds_RaList_h
