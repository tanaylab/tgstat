#ifndef STATQUADTREECACHEDSERIALIZER_H_
#define STATQUADTREECACHEDSERIALIZER_H_

#include <math.h>
#include <strings.h>

#include "StatQuadTreeCached.h"

// StatQuadTreeCachedSerializer is an auxiliary class that helps to construct and serialize StatQuadTreeCached.
// StatQuadTreeCached can be constructed by two different ways:
// 1. By constructing first StatQuadTree and then passing it to StatQuadTreeCached::serialize.
// 2. By using StatQuadTreeCachedSerializer.
// 
// The first ("standard") way is used for the majority of cases. However when the number of objects is too high to hold them all
// in memory (they must be first inserted into StatQuadTree), StatQuadTreeCachedSerializer may come handy.
// 
// StatQuadTreeCachedSerializer allows construction of enormously large StatQuadTreeCached by splitting the inserting the objects
// into a designated "subtrees". Here's how it works:
// 
// 1. Call begin() and specify the arena and the number of subtrees. This number must be a power of 4. For evenly distrubuted objects
//    over the arena the total memory consumption will be reduced by the number of subtrees.
// 2. begin() splits the arena to subarenas (can be retrieved by get_subarenas()). Each object lays inside one or more subarenas
//    (object can cross the border between subarenas and thus belong to a few subarenas).
// 3. Objects can now be inserted using insert() function. Yet an important condition must be met:
//    OBJECTS INSERTION MUST BE SORTED BY SUBARENAS!
//    It means that all the objects belonging to a certain subarena must be added one after another without being interrupted by objects
//    belonging to a different subarena. Once an object belonging to another subarena is inserted the last subarena is sealed and no more objects
//    can be inserted to it. Violation of this condition will generate an error.
// 4. Object insertion is finished by calling end() function.
//
// An object inserted to StatQuadTreeCachedSerializer must have the following functions:
//    int64_t get_x1() const;
//    int64_t get_y1() const;
//    int64_t get_x2() const;
//    int64_t get_y2() const;

template <class T, class Size>
class StatQuadTreeCachedSerializer {
public:
	typedef T ValueType;

	StatQuadTreeCachedSerializer();

	~StatQuadTreeCachedSerializer() { end(); }

	// Note: num_subtrees must be a power of 4
	void begin(BufferedFile &file, int64_t x1, int64_t y1, int64_t x2, int64_t y2, unsigned num_subtrees, bool check_overlaps, 
			   int64_t chunk_size, int64_t max_num_chunks, unsigned max_depth = 20, unsigned max_node_objs = 20);

	void end();

	// Note: OBJECTS INSERTION MUST BE SORTED BY SUBARENAS
	void insert(const T &obj);
	
	const Rectangles &get_subarenas() const { return m_subarenas; }

	// returns the index of the last subtree used for object insertion (can be -1 if the objects fall at the border between the subtrees)
	int get_cur_subarena_idx() const { return m_cur_qtree_idx; }

private:
	BufferedFile            *m_file;
	bool                     m_is_ended;
	Rectangle                m_arena;
	unsigned                 m_num_subtrees;
	unsigned                 m_num_subtrees_sqrt;
	bool                     m_check_overlaps;
	int64_t                  m_chunk_size;
	int64_t                  m_max_num_chunks;
	unsigned                 m_max_depth;
	unsigned                 m_max_node_objs;
	int64_t                  m_tree_start_fpos;
	int64_t                  m_top_chunk_start_fpos;
	Rectangles               m_subarenas;
	vector<int64_t>          m_subtree_fpos;
	vector<bool>             m_is_qtree_sealed;
	vector<typename StatQuadTreeCached<T, Size>::Stat> m_stat;

	Size                     m_cur_id;
	Size                     m_cur_id_offset;   // equals to number of unique objects in the already seal quad trees
	StatQuadTree<T, Size>    m_cur_qtree;
	int                      m_cur_qtree_idx;
	vector<T>                m_border_objs;
	vector<Size>             m_border_obj_ids;
	vector< vector<size_t> > m_border_obj_ptrs;

	int idx2dto1d(int i, int j) const { return i + m_num_subtrees_sqrt * j; }
	Rectangle &subarena(int i, int j) { return m_subarenas[idx2dto1d(i, j)]; }

	void set_subarenas(int i1, int j1, int i2, int j2, int64_t x1, int64_t y1, int64_t x2, int64_t y2);

	// serializes the top tree (the ancestor of all subtrees); returns file position relative to the beginning of the chunk
	int64_t serialize_top_tree(int i1, int j1, int i2, int j2, int64_t x1, int64_t y1, int64_t x2, int64_t y2,
							   typename StatQuadTreeCached<T, Size>::Stat &stat);

	void seal_qtree();
};

typedef StatQuadTreeCachedSerializer<Rectangle_val<float>, uint64_t> RectsQuadTreeCachedSerializer;
typedef StatQuadTreeCachedSerializer<Point_val<float>, uint64_t>     PointsQuadTreeCachedSerializer;

//------------------------------ IMPLEMENTATION ----------------------------------------

template <class T, class Size>
StatQuadTreeCachedSerializer<T, Size>::StatQuadTreeCachedSerializer() :
	m_file(NULL), m_is_ended(true), m_num_subtrees(0), m_num_subtrees_sqrt(0), m_cur_qtree_idx(-1)
{}

template <class T, class Size>
void StatQuadTreeCachedSerializer<T, Size>::begin(BufferedFile &file, int64_t x1, int64_t y1, int64_t x2, int64_t y2, unsigned num_subtrees, bool check_overlaps, 
												  int64_t chunk_size, int64_t max_num_chunks, unsigned max_depth, unsigned max_node_objs)
{
	m_file = &file;
	m_cur_id = 0;
	m_cur_id_offset = 0;
	m_cur_qtree_idx = -1;
	m_check_overlaps = check_overlaps;
	m_chunk_size = chunk_size;
	m_max_num_chunks = max_num_chunks;
	m_max_depth = max_depth;
	m_max_node_objs = max_node_objs;
	m_border_objs.clear();
	m_border_obj_ids.clear();

	m_arena.x1 = x1;
	m_arena.y1 = y1;
	m_arena.x2 = x2;
	m_arena.y2 = y2;

	// check that num_subtrees is a power of 4:
	// the position of the most significant must be odd + the rest of the bits must be zero
	int msbit_pos = ffs(num_subtrees);
	if (num_subtrees != 1 << (msbit_pos - 1) || msbit_pos % 2 == 0)
		TGLError< StatQuadTreeCachedSerializer<T, Size> >("Number of sub quad trees must be a power of 4");

	m_num_subtrees = num_subtrees;
	m_num_subtrees_sqrt = (unsigned)(sqrt(m_num_subtrees) + .5);

	m_subarenas.resize(m_num_subtrees);
	m_subtree_fpos.resize(m_num_subtrees);
	m_is_qtree_sealed.resize(m_num_subtrees, false);
	m_border_obj_ptrs.resize(m_num_subtrees);
	m_stat.resize(m_num_subtrees);

	set_subarenas(0, 0, m_num_subtrees_sqrt, m_num_subtrees_sqrt, m_arena.x1, m_arena.y1, m_arena.x2, m_arena.y2);

	if (m_num_subtrees > 1) {
		m_tree_start_fpos = m_file->tell();

		uint64_t num_objs = 0;
		int64_t root_chunk_start_fpos = 0;
		m_file->write(&num_objs, sizeof(num_objs));  // write a placeholder for number of objects
		m_file->write(&root_chunk_start_fpos, sizeof(root_chunk_start_fpos)); // write a placeholder for the file position of the root chunk
	}

	m_is_ended = false;
}

template <class T, class Size>
void StatQuadTreeCachedSerializer<T, Size>::end()
{
	if (m_is_ended) 
		return;

	m_is_ended = true;

	if (m_num_subtrees == 1) {
		StatQuadTreeCached<T, Size> qtree_cached(m_chunk_size, m_max_num_chunks);
		qtree_cached.serialize(*m_file, m_cur_qtree);
	} else {
		seal_qtree();

		// seal ALL unsealed tree (there might be some empty ones or those which objects are at the border)
		for (unsigned i = 0; i < m_num_subtrees; ++i) {
			if (!m_is_qtree_sealed[i]) {
				m_cur_qtree_idx = i;
				m_cur_qtree.init(m_subarenas[i].x1, m_subarenas[i].y1, m_subarenas[i].x2, m_subarenas[i].y2, m_max_depth, m_max_node_objs);
				seal_qtree();
			}
		}

		uint64_t num_objs = (uint64_t)m_cur_id_offset;
		int64_t root_chunk_start_fpos = m_file->tell();

		// write down the header of the tree
		m_file->seek(m_tree_start_fpos, SEEK_SET);
		m_file->write(&num_objs, sizeof(num_objs));
		m_file->write(&root_chunk_start_fpos, sizeof(root_chunk_start_fpos));
		m_file->seek(root_chunk_start_fpos, SEEK_SET);

		if (num_objs) {
			int64_t top_node_offset = 0;

			// The header of the chunk starts with the chunk size in bytes. To calculate it we need to know the number of nodes in the chunk.
			// It equals to 1 + 4 + 4^2 + ... + 4^n, where n is the depth of the tree. This equals to (1 - 4^n)/(1-4) (= sum of geometrical progression).
			// The number of subtrees (= m_num_subtrees) is known and it is equal to 4^n.
			// Therefore the total number of nodes equals to (m_num_subtrees - 1) / 3.
			int64_t size = sizeof(typename StatQuadTreeCached<T, Size>::Node) * (m_num_subtrees - 1) / 3 + sizeof(size) + sizeof(top_node_offset);

			m_top_chunk_start_fpos = m_file->tell();
			m_file->write(&size, sizeof(size));                       // write the chunk size
			m_file->write(&top_node_offset, sizeof(top_node_offset)); // write a placeholder for a top node file position

			typename StatQuadTreeCached<T, Size>::Stat stat;
			top_node_offset = serialize_top_tree(0, 0, m_num_subtrees_sqrt, m_num_subtrees_sqrt, m_arena.x1, m_arena.y1, m_arena.x2, m_arena.y2, stat);

			int64_t cur_fpos = m_file->tell();
			m_file->seek(m_top_chunk_start_fpos + sizeof(size), SEEK_SET);
			m_file->write(&top_node_offset, sizeof(top_node_offset));
			m_file->seek(cur_fpos, SEEK_SET);
		}
	}
}

template <class T, class Size>
void StatQuadTreeCachedSerializer<T, Size>::insert(const T &obj)
{
	// if object does not belong to the current quad tree
	if (m_cur_qtree_idx < 0 || !obj.do_intersect(m_subarenas[m_cur_qtree_idx])) {
		unsigned new_qtree_idx;

		// find to which one it belongs to
		for (new_qtree_idx = 0; new_qtree_idx < m_num_subtrees; ++new_qtree_idx) {
			if (obj.is_inside(m_subarenas[new_qtree_idx]))
				break;
		}

		// seal the old qtree
		if (new_qtree_idx < m_num_subtrees) {
			if (m_cur_qtree_idx >= 0)
				seal_qtree();
			m_cur_qtree_idx = new_qtree_idx;
			m_cur_qtree.init(m_subarenas[m_cur_qtree_idx].x1, m_subarenas[m_cur_qtree_idx].y1, m_subarenas[m_cur_qtree_idx].x2, m_subarenas[m_cur_qtree_idx].y2,
							 m_max_depth, m_max_node_objs);
		}
	}

	if (m_cur_qtree_idx >= 0 && obj.is_inside(m_subarenas[m_cur_qtree_idx])) { // object is inside the quad tree
		if (m_is_qtree_sealed[m_cur_qtree_idx])
			TGLError< StatQuadTreeCachedSerializer<T, Size> >("Objects are inserted to StatQuadTreeCachedSerializer unordered");

		if (m_check_overlaps && m_cur_qtree.do_intersect(Rectangle(obj.get_x1(), obj.get_y1(), obj.get_x2(), obj.get_y2())))
			TGLError< StatQuadTreeCachedSerializer<T, Size> >("Inserted object (%ld, %ld)-(%ld, %ld) intersects existing ones",
															  obj.get_x1(), obj.get_y1(), obj.get_x2(), obj.get_y2());

		m_cur_qtree.insert(obj);
		m_cur_id++;
	} else {
		bool inserted = false;

		for (unsigned i = 0; i < m_num_subtrees; ++i) {
			if (obj.do_intersect(m_subarenas[i])) {
				if (m_is_qtree_sealed[i])
					TGLError< StatQuadTreeCachedSerializer<T, Size> >("Objects are inserted to StatQuadTreeCachedSerializer unordered");

				if (!inserted) {
					m_border_objs.push_back(obj);
					m_border_obj_ids.push_back(-1);
					inserted = true;
				}
				m_border_obj_ptrs[i].push_back(m_border_objs.size() - 1);
			}
		}
	}
}

template <class T, class Size>
void StatQuadTreeCachedSerializer<T, Size>::seal_qtree()
{
	if (m_cur_qtree_idx < 0 || m_num_subtrees <= 1)
		return;

	// insert all the objects that are laying at the border
	uint64_t num_unique_objs = m_cur_qtree.get_num_objs();
	vector<Size> local2global_id(num_unique_objs);

	for (uint64_t i = 0; i < num_unique_objs; ++i) 
		local2global_id[i] = i + m_cur_id_offset;

	for (vector<size_t>::const_iterator iptr = m_border_obj_ptrs[m_cur_qtree_idx].begin(); iptr != m_border_obj_ptrs[m_cur_qtree_idx].end(); ++iptr) {
		if (m_border_obj_ids[*iptr] == -1) { // object hasn't been added yet to any quad tree
			m_border_obj_ids[*iptr] = m_cur_id_offset + num_unique_objs;
			++num_unique_objs;
		}

		const T &obj = m_border_objs[*iptr];

		if (m_check_overlaps && m_cur_qtree.do_intersect(Rectangle(obj.get_x1(), obj.get_y1(), obj.get_x2(), obj.get_y2())))
			TGLError< StatQuadTreeCachedSerializer<T, Size> >("Inserted object (%ld, %ld)-(%ld, %ld) intersects existing ones",
															  obj.get_x1(), obj.get_y1(), obj.get_x2(), obj.get_y2());

		m_cur_qtree.insert(obj);
		local2global_id.push_back(m_border_obj_ids[*iptr]);
	}

	StatQuadTreeCached<T, Size> qtree_cached(m_chunk_size, m_max_num_chunks);
	m_subtree_fpos[m_cur_qtree_idx] = qtree_cached.serialize_as_chunk(*m_file, m_cur_qtree, local2global_id);

	m_cur_id_offset += num_unique_objs;
	m_stat[m_cur_qtree_idx] = m_cur_qtree.m_nodes.front().stat;
	m_cur_qtree.reset();
	m_is_qtree_sealed[m_cur_qtree_idx] = true;
}

template <class T, class Size>
void StatQuadTreeCachedSerializer<T, Size>::set_subarenas(int i1, int j1, int i2, int j2, int64_t x1, int64_t y1, int64_t x2, int64_t y2)
{
	if (x1 == x2 || y1 == y2)
		TGLError< StatQuadTreeCachedSerializer<T, Size> >("Arena is not big enough to be split to %ld subtrees", m_subarenas.size());
	
	if (i1 >= i2 - 1) {
		Rectangle &rect = subarena(i1, j1);
		rect.x1 = x1;
		rect.y1 = y1;
		rect.x2 = x2;
		rect.y2 = y2;
	} else {
		int64_t split_x = (x1 + x2) / 2;
		int64_t split_y = (y1 + y2) / 2;
		int split_i = (i1 + i2) / 2;
		int split_j = (j1 + j2) / 2;

		set_subarenas(i1, j1, split_i, split_j, x1, y1, split_x, split_y);  // SW
		set_subarenas(split_i, j1, i2, split_j, split_x, y1, x2, split_y);  // SE
		set_subarenas(i1, split_j, split_i, j2, x1, split_y, split_x, y2);  // NW
		set_subarenas(split_i, split_j, i2, j2, split_x, split_y, x2, y2);  // NE
	}
}

template <class T, class Size>
int64_t StatQuadTreeCachedSerializer<T, Size>::serialize_top_tree(int i1, int j1, int i2, int j2, int64_t x1, int64_t y1, int64_t x2, int64_t y2, typename StatQuadTreeCached<T, Size>::Stat &stat)
{
	typename StatQuadTreeCached<T, Size>::Node node;

	node.is_leaf = false;
	node.arena = Rectangle(x1, y1, x2, y2);

	if (i2 - i1 <= 2) {  // lowest node in the top tree (== below are only serialized already subtrees)
		for (int iquad = 0; iquad < StatQuadTreeCached<T, Size>::NUM_QUADS; iquad++) {
			int idx;

			if (iquad == StatQuadTreeCached<T, Size>::SW)
				idx = idx2dto1d(i1, j1);
			else if (iquad == StatQuadTreeCached<T, Size>::SE)
				idx = idx2dto1d(i1 + 1, j1);
			else if (iquad == StatQuadTreeCached<T, Size>::NW)
				idx = idx2dto1d(i1, j1 + 1);
			else if (iquad == StatQuadTreeCached<T, Size>::NE)
				idx = idx2dto1d(i1 + 1, j1 + 1);

			node.kid_ptr[iquad] = -m_subtree_fpos[idx]; // lowest node in the top tree points to the absolute position of chunk in a file (marked by negative number)
			node.stat.weighted_sum += m_stat[idx].weighted_sum;
			node.stat.min_val = min(node.stat.min_val, m_stat[idx].min_val);
			node.stat.max_val = max(node.stat.max_val, m_stat[idx].max_val);
			node.stat.occupied_area += m_stat[idx].occupied_area;
		}
	} else {
		int64_t split_x = (x1 + x2) / 2;
		int64_t split_y = (y1 + y2) / 2;
		int split_i = (i1 + i2) / 2;
		int split_j = (j1 + j2) / 2;

		for (int iquad = 0; iquad < StatQuadTreeCached<T, Size>::NUM_QUADS; iquad++) {
			if (iquad == StatQuadTreeCached<T, Size>::SW) 
				node.kid_ptr[iquad] = serialize_top_tree(i1, j1, split_i, split_j, x1, y1, split_x, split_y, node.stat);  // SW
			else if (iquad == StatQuadTreeCached<T, Size>::SE)
				node.kid_ptr[iquad] = serialize_top_tree(split_i, j1, i2, split_j, split_x, y1, x2, split_y, node.stat);  // SE
			else if (iquad == StatQuadTreeCached<T, Size>::NW)
				node.kid_ptr[iquad] = serialize_top_tree(i1, split_j, split_i, j2, x1, split_y, split_x, y2, node.stat);  // NW
			else if (iquad == StatQuadTreeCached<T, Size>::NE)
				node.kid_ptr[iquad] = serialize_top_tree(split_i, split_j, i2, j2, split_x, split_y, x2, y2, node.stat);  // NE
		}
	}

	int64_t top_node_fpos = m_file->tell();
	m_file->write(&node, sizeof(node));
	stat.weighted_sum += node.stat.weighted_sum;
	stat.min_val = min(stat.min_val, node.stat.min_val);
	stat.max_val = max(stat.max_val, node.stat.max_val);
	stat.occupied_area += node.stat.occupied_area;

	return top_node_fpos - m_top_chunk_start_fpos;
}

#endif

