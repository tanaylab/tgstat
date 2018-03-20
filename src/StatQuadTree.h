/*
 * StatQuadTree.h
 *
 *  Created on: Dec 19, 2011
 *      Author: hoichman
 */

#ifndef QUADTREE_H_
#define QUADTREE_H_

#include <stdio.h>
#include <algorithm>
#include <limits>
#include <math.h>
#include <queue>
#include <vector>

#include "BufferedFile.h"
#include "DiagonalBand.h"
#include "Point.h"
#include "Rectangle.h"

using namespace std;

#pragma pack(push)
#pragma pack(8)

//----------------------------------------- Rectangle_val -----------------------------------

template<typename T>
struct Rectangle_val : public Rectangle {
	T v;

	Rectangle_val() {}
	Rectangle_val(int64_t _x1, int64_t _y1, int64_t _x2, int64_t _y2, const T &_v = T()) : Rectangle(_x1, _y1, _x2, _y2), v(_v) {}
	Rectangle_val(const Rectangle &rect, const T &_v = T()) : Rectangle(rect), v(_v) {}

	double val(const Rectangle &, void *) const { return v; }
	double val(const Rectangle &, const DiagonalBand &, void *) const { return v; }

	char *debug_str() const {
		static char str[200];
		sprintf(str, "(%lld - %lld) (%lld - %lld) %g", x1, x2, y1, y2, (double)v);
		return str;
	}
};

//----------------------------------------- Point_val ---------------------------------------

template<typename T>
struct Point_val : public Point {
	T v;

	Point_val() {}
	Point_val(int64_t _x, int64_t _y, const T &_v = T()) : Point(_x, _y), v(_v) {}
	Point_val(const Point &point, const T &_v = T()) : Point(point), v(_v) {}

	double val(const Rectangle &, void *) const { return v; }
	double val(const Rectangle &, const DiagonalBand &, void *) const { return v; }

	char *debug_str() const {
		static char str[200];
		sprintf(str, "(%ld - %ld) %g", x, y, (double)v);
		return str;
	}
};

#pragma pack(pop)

// StatQuadTree stores non-overlapping geographical objects and their value and allows efficient retrieval of various
// statistics about these objects.
//
// StatQuadTree must be initialized with coordinates of a "master-rectangle" ("arena") that contains all the inserted objects AND
// the area of the queries!
//
// !!!! NOTE: The objects stored in StatQuadTree (i.e. "T") must be disjoint, i.e. not to intersect with each other!!!!
//            An attempt to insert overlapping objects might result in degenerative data structure and excessive memory usage.
//            (Overlapping objects might trigger never-ending quad splitting that is bound only by max_depth.)
// 
// "Size" template argument indicates the type that is used to store object ids. Generally it is expected to be either unsigned int
// when the number of objects is lower than 2^32 or uint64_t if the number of objects exceeds 2^32.
//
// An object stored in StatQuadTree must have the following functions:
//    // returns true if the object intersects with rect
//    bool do_intersect(const Rectangle &rect) const;
//
//    // returns the value associated with the object for the given area
//    double val(const Rectangle &area, void *user_ptr) const;
//
// For using intersect() function the object must also have:
//    // returns an area (rectangle) of an intersection between obj and rect
//    Rectangle intersect(const Rectangle &rect) const;
//
// For using debug_print_tree() function the object must also have:
//    // returns a string that represens the object for debug print outs
//    char *debug_str() const;
//
// For using serialize() / unserialize() functions the object must support primitive serialization: i.e.:
//    memcpy(dst, &obj, sizeof(obj))
//  It means that the object cannot contain pointers or any other objects that contain pointers inside (like STL containers).

template <class T, class Size>
class StatQuadTree {
public:
	template<class U, class V> friend class StatQuadTreeCached;
	template<class U, class V> friend class StatQuadTreeCachedSerializer;

	typedef T ValueType;

#pragma pack(push)
#pragma pack(8)

	struct Stat {
		int64_t occupied_area;
		double  weighted_sum;
		double  min_val;
		double  max_val;

		Stat() : occupied_area(0), weighted_sum(0.), min_val(numeric_limits<double>::max()), max_val(-numeric_limits<double>::max()) {}

		void reset() {
			occupied_area = 0;
			weighted_sum = 0.;
			min_val = numeric_limits<double>::max();
			max_val = -numeric_limits<double>::max();
		}
	};

#pragma pack(pop)

	StatQuadTree() : m_uptr(NULL) { init(0, 0, 0, 0); }

	// x1, y1, x2, y2 - are the coordinates of the "arena" - the rectangle that should contain all the inserted objects.
	// max_depth      - maximal quad-tree depth. The maximal number of nodes would be 4^max_depth.
	// max_node_objs  - maximal number of objects in the node. Once the number of objects exceed max_node_objs, the node is split.
	StatQuadTree(int64_t x1, int64_t y1, int64_t x2, int64_t y2, unsigned max_depth = 20, unsigned max_node_objs = 20);

	void init(int64_t x1, int64_t y1, int64_t x2, int64_t y2, unsigned max_depth = 20, unsigned max_node_objs = 20);

	void insert(const T &obj);

	// returns true if rect intersects any of the rectangles in the quad-tree
	bool do_intersect(const Rectangle &rect) const;

	// intersects rect with the objects stored in the quad-tree, returns a vector of intersected areas and
	// the indices of the objects that created these intersected areas
	void intersect(const Rectangle &rect, Rectangles &intersection, vector<Size> &intersected_objs_indices) const;

	// same as a bandless intersect() version, however the returned rectangles are shrinked to the band (see DiagonalBand::shrink2intersected)
	void intersect(const Rectangle &rect, const DiagonalBand &band, Rectangles &intersection, vector<Size> &intersected_objs_indices) const;

	// returns statitstics for the given rectangle: weighted sum, occupied area, minimal and maximal value of the object that intersect with rect.
	// Note: the average value can be achieved by dividing weighted sum by occupied.area
	void get_stat(const Rectangle &rect, Stat &stat) const;

	// returns a "density map" of the quad tree, i.e. all the rectangles covered by the leaves of the quad tree
	void get_density_map(Rectangles &density_map) const;

	void debug_print_tree() const;

	void reset();

	void reset(int64_t x1, int64_t y1, int64_t x2, int64_t y2);

	// prior to serialize() file must be open for writing
	void serialize(BufferedFile &file);

	// prior to unserialize() file must be open for reading
	void unserialize(BufferedFile &file);

	uint64_t get_num_objs() const { return m_objs.size(); }

	bool empty() const { return m_objs.empty(); }

	const Rectangle &get_arena() const { return m_nodes.front().arena; }

	unsigned get_max_depth() const { return m_max_depth; }

	unsigned get_max_node_objs() const { return m_max_node_objs; }

	const vector<T> &get_objs() const { return m_objs; }

	void set_uptr(void *uptr) { m_uptr = uptr; }

private:
	// nw = North West, ne = North East, se = South East, sw = South West;
	enum { NW, NE, SE, SW, NUM_QUADS };

#pragma pack(push)
#pragma pack(8)

	struct Node {
		struct Node_data {
			Size kid_idx[NUM_QUADS];
		};

		struct Leaf_data {
			// since each object might be split between several leaves, the number of total sub-objects might exceed 2^32 => use uint64_t
			uint64_t obj_ptr_start_idx;
			uint64_t obj_ptr_end_idx;
		};

		union {
			Node_data node;
			Leaf_data leaf;
		};

		bool   is_leaf;
		Stat   stat;
		Rectangle arena;

		Node() {}

		Node(const Rectangle &_arena) : arena(_arena) {
			init_leaf();
		}

		void init_leaf() {
			is_leaf = true;
			leaf.obj_ptr_start_idx = leaf.obj_ptr_end_idx = 0;
		}

		void init_node() {
			is_leaf = false;
			node.kid_idx[NW] = node.kid_idx[NE] = node.kid_idx[SE] = node.kid_idx[SW] = -1;
		}
	};

#pragma pack(pop)

	vector<Node>          m_nodes;
	vector<Size>          m_obj_ptrs;
	vector<uint64_t>      m_obj_ptrs_free_chunks;
	vector<T>             m_objs;
	mutable vector<bool>  m_intersected_objs; // used for intersection queries, not serialized
	unsigned              m_max_depth;
	unsigned              m_max_node_objs;
	void                 *m_uptr;

	void insert(Node *&node, const Rectangle &intersection, unsigned depth, const T &obj, uint64_t obj_idx);
	void insert2leaf(Node *&node, uint64_t obj_idx);
	void create_quad(Node *&node, int quad, const Rectangle &arena);
	bool do_intersect(const Node &node, const Rectangle &rect) const;
	void intersect(const Node &node, const Rectangle &rect, Rectangles &intersection, vector<Size> &intersected_objs_indices) const;
	void intersect(const Node &node, const Rectangle &rect, const DiagonalBand &band, Rectangles &intersection, vector<Size> &intersected_objs_indices) const;
	void get_stat(const Node &node, const Rectangle &rect, Stat &stat) const;
	void update_stat(const T &obj, Stat &stat, const Rectangle &intersection) const;
	void update_stat(const Node &node, Stat &stat) const;
	void debug_print_tree(const Node &node, unsigned depth) const;

public:
	// Nearest neighbor iterator returns the objects according to their proximity (Manhattan distance)
	// to query rectangle starting from the nearest one.
	// Use excluded area to exclude objects that intersect that area.
	// Unlike standard STL iterators you can't copy this iterator nor compare it to another one.
	// Call is_end() to test whether the iterator reached the end.
	class NNIterator {
	public:
		NNIterator(StatQuadTree<T, Size> *parent) : m_parent(parent) {}

		bool begin(const Rectangle &query, Rectangle excluded_area = Rectangle(0, 0, -1, -1)); // returns true if end is not reached yet
		bool next();                                                                           // returns true if end is not reached yet
		bool is_end() const { return m_neighbors.empty(); }

		T &operator*() { return cur_obj(); }
		const T &operator*() const { return cur_obj(); }

		T *operator->() { return &cur_obj(); }
		const T *operator->() const { return &cur_obj(); }

	private:
		struct Neighbor {
			Node    *node;
			T       *obj;
			int64_t  dist;

			Neighbor(Node *_node, T *_obj, int64_t _dist) : node(_node), obj(_obj), dist(_dist) {}
			Neighbor(const Neighbor &_obj) : node(_obj.node), obj(_obj.obj), dist(_obj.dist) {}

			// The order in the heap must give preference to objects before nodes.
			// Otherwise nearest neighbor query that covers the whole arena will cause all the tree nodes to be added to the heap.
			bool operator<(const Neighbor &obj) const { return dist > obj.dist || (dist == obj.dist && node); }
		};

		Rectangle                 m_query;
		Rectangle                 m_excluded_area;
		StatQuadTree<T, Size>    *m_parent;
		priority_queue<Neighbor>  m_neighbors;
		vector<bool>              m_used_objs;  // used to prevent multiple return of the same object that belongs to different quads

		T &cur_obj() { return *m_neighbors.top().obj; }
		const T &cur_obj() const { return *m_neighbors.top().obj; }

		// block copy constructor and default operators
		NNIterator(const NNIterator &) {}
		NNIterator &operator=(const NNIterator &) { return *this; }
		bool operator==(const NNIterator &) { return false; }
	};
};


typedef StatQuadTree<Rectangle_val<float>, uint64_t> RectsQuadTree;
typedef StatQuadTree<Point_val<float>, uint64_t>     PointsQuadTree;

//------------------------------ IMPLEMENTATION ----------------------------------------

//============================== NNIterator =============================================

template <class T, class Size>
bool StatQuadTree<T, Size>::NNIterator::begin(const Rectangle &query, Rectangle excluded_area)
{
	m_query = query;
	m_excluded_area = excluded_area;
	m_used_objs.clear();
	m_used_objs.resize(m_parent->m_objs.size(), false);
	m_neighbors = priority_queue<Neighbor>();

	if (!m_parent->m_nodes.size())
		return false;

	Node &root = m_parent->m_nodes.front();
	if (!root.arena.is_inside(m_excluded_area)) 
		m_neighbors.push(Neighbor(&root, NULL, root.arena.manhattan_dist(m_query, true)));
	return next();
}

template <class T, class Size>
bool StatQuadTree<T, Size>::NNIterator::next()
{
	if (m_neighbors.empty())
		return false;

	if (m_neighbors.top().obj)
		m_neighbors.pop();

	while (!m_neighbors.empty()) {
		if (m_neighbors.top().obj)
			return true;

		Node *node = m_neighbors.top().node;
		m_neighbors.pop();
		if (node->is_leaf) {
			for (uint64_t iobj_ptr = node->leaf.obj_ptr_start_idx; iobj_ptr < node->leaf.obj_ptr_end_idx; ++iobj_ptr) {
				size_t idx = m_parent->m_obj_ptrs[iobj_ptr];
				if (!m_used_objs[idx]) {
					T &obj = m_parent->m_objs[idx];
					if (!obj.do_intersect(m_excluded_area)) {
						m_neighbors.push(Neighbor(NULL, &obj, obj.manhattan_dist(m_query, true)));
						m_used_objs[idx] = true;
					}
				}
			}
		} else {
			for (int iquad = 0; iquad < NUM_QUADS; ++iquad) {
				Node *quad = &m_parent->m_nodes[node->node.kid_idx[iquad]];
				if (quad->stat.occupied_area > 0 && !quad->arena.is_inside(m_excluded_area))
					m_neighbors.push(Neighbor(quad, NULL, quad->arena.manhattan_dist(m_query, true)));
			}
		}
	}
	return false;
}


//============================== StatQuadTree ===========================================

template <class T, class Size>
StatQuadTree<T, Size>::StatQuadTree(int64_t x1, int64_t y1, int64_t x2, int64_t y2, unsigned max_depth, unsigned max_node_objs) :
	m_uptr(NULL)
{
	init(x1, y1, x2, y2, max_depth, max_node_objs);
}

template <class T, class Size>
void StatQuadTree<T, Size>::init(int64_t x1, int64_t y1, int64_t x2, int64_t y2, unsigned max_depth, unsigned max_node_objs)
{
	m_max_depth = max_depth;
	m_max_node_objs = max_node_objs;

	reset(x1, y1, x2, y2);
}

template <class T, class Size>
void StatQuadTree<T, Size>::reset()
{
	Node &root = m_nodes.front();
	reset(root.arena.x1, root.arena.y1, root.arena.x2, root.arena.y2);
}

template <class T, class Size>
void StatQuadTree<T, Size>::reset(int64_t x1, int64_t y1, int64_t x2, int64_t y2)
{
	m_nodes.clear();
	m_obj_ptrs.clear();
	m_objs.clear();
	m_obj_ptrs_free_chunks.clear();
	m_intersected_objs.clear();

	// allocate the root
	m_nodes.push_back(Node(Rectangle(x1, y1, x2, y2)));
}

template <class T, class Size>
void StatQuadTree<T, Size>::insert(const T &obj)
{
	m_objs.push_back(obj);
	Rectangle intersection = obj.intersect(m_nodes.front().arena);
	if (intersection.is_non_empty_area()) {
		Node *proot = &m_nodes.front();
		insert(proot, intersection, 0, obj, m_objs.size() - 1);
	}
}

template <class T, class Size>
void StatQuadTree<T, Size>::insert(Node *&node, const Rectangle &intersection, unsigned depth, const T &obj, uint64_t obj_idx)
{
	// update node statistics
	update_stat(obj, node->stat, intersection);

	// add the object itself
	if (node->is_leaf) {
		if (node->leaf.obj_ptr_end_idx - node->leaf.obj_ptr_start_idx < m_max_node_objs || depth >= m_max_depth ||
				node->arena.width() < 4 || node->arena.height() < 4)
		{
			insert2leaf(node, obj_idx);
			return;
		}

		uint64_t obj_ptr_start_idx = node->leaf.obj_ptr_start_idx;
		uint64_t obj_ptr_end_idx = node->leaf.obj_ptr_end_idx;

		// convert the leaf to node
		node->init_node();

		int64_t split_x = (node->arena.x1 + node->arena.x2) / 2;
		int64_t split_y = (node->arena.y1 + node->arena.y2) / 2;

		create_quad(node, NW, Rectangle(node->arena.x1, split_y, split_x, node->arena.y2));
		create_quad(node, NE, Rectangle(split_x, split_y, node->arena.x2, node->arena.y2));
		create_quad(node, SE, Rectangle(split_x, node->arena.y1, node->arena.x2, split_y));
		create_quad(node, SW, Rectangle(node->arena.x1, node->arena.y1, split_x, split_y));

		// reassign all objects to newly created quads
		for (uint64_t i = obj_ptr_start_idx; i < obj_ptr_end_idx; i++) {
			for (int iquad = 0; iquad < NUM_QUADS; iquad++) {
				Node *quad = &m_nodes[node->node.kid_idx[iquad]];
				Rectangle intersection = m_objs[m_obj_ptrs[i]].intersect(quad->arena);
				if (intersection.is_non_empty_area())
					insert(quad, intersection, depth + 1, m_objs[m_obj_ptrs[i]], m_obj_ptrs[i]);
			}
		}

		m_obj_ptrs_free_chunks.push_back(obj_ptr_start_idx);
	}

	// add an object to the quads
	for (int iquad = 0; iquad < NUM_QUADS; iquad++) {
		Node *quad = &m_nodes[node->node.kid_idx[iquad]];
		Rectangle intersection = obj.intersect(quad->arena);
		if (intersection.is_non_empty_area()) {
			// after insert() the node pointer might change (m_nodes might get resized)
			size_t node_idx = node - &*m_nodes.begin();
			insert(quad, intersection, depth + 1, obj, obj_idx);
			node = &m_nodes[node_idx];
		}
	}
}

template <class T, class Size>
void StatQuadTree<T, Size>::insert2leaf(Node *&node, uint64_t obj_idx)
{
	// allocate a chunk in m_obj_ptrs if needed
	if (node->leaf.obj_ptr_end_idx == node->leaf.obj_ptr_start_idx) {
		if (m_obj_ptrs_free_chunks.empty()) {
			uint64_t old_size = m_obj_ptrs.size();
			m_obj_ptrs.resize(old_size + m_max_node_objs);
			node->leaf.obj_ptr_start_idx = node->leaf.obj_ptr_end_idx = old_size;
		} else {
			node->leaf.obj_ptr_start_idx = node->leaf.obj_ptr_end_idx = m_obj_ptrs_free_chunks.back();
			m_obj_ptrs_free_chunks.pop_back();
		}
	} else if (node->leaf.obj_ptr_end_idx - node->leaf.obj_ptr_start_idx >= m_max_node_objs) {
		uint64_t num_objs = node->leaf.obj_ptr_end_idx - node->leaf.obj_ptr_start_idx;
		unsigned num_chunks = num_objs / m_max_node_objs;

		// We are going to add the object to the leaf anyway even if the number of objects exceeds m_max_node_objs.
		// Check that we have enough space for the insert. Our space is always a power of 2 multiplied by m_max_node_objs, i.e.:
		// m_max_node_objs * 2^n.
		//
		// ...So here's the test. It's a bit tricky and uses ffs() that returns the position of the most significant bit.
		if (m_max_node_objs << (ffs(num_chunks) - 1) == num_objs) {
			// OK, there's no space for one more element, we have to reallocate the space
			uint64_t old_size = m_obj_ptrs.size();

			m_obj_ptrs.resize(old_size + 2 * num_objs);
			copy(m_obj_ptrs.begin() + node->leaf.obj_ptr_start_idx, m_obj_ptrs.begin() + node->leaf.obj_ptr_end_idx, m_obj_ptrs.begin() + old_size);

			for (unsigned i = 0; i < num_chunks; i++)
				m_obj_ptrs_free_chunks.push_back(node->leaf.obj_ptr_start_idx + i * m_max_node_objs);

			node->leaf.obj_ptr_start_idx = old_size;
			node->leaf.obj_ptr_end_idx = old_size + num_objs;
		}
	}

	m_obj_ptrs[node->leaf.obj_ptr_end_idx] = obj_idx;
	node->leaf.obj_ptr_end_idx++;
}

template <class T, class Size>
void StatQuadTree<T, Size>::create_quad(Node *&node, int quad, const Rectangle &arena)
{
	size_t node_idx = node - &*m_nodes.begin();
	node->node.kid_idx[quad] = m_nodes.size();
	m_nodes.push_back(Node(arena));
	node = &m_nodes[node_idx];
}

template <class T, class Size>
bool StatQuadTree<T, Size>::do_intersect(const Rectangle &rect) const
{
	return do_intersect(m_nodes.front(), rect);
}

template <class T, class Size>
bool StatQuadTree<T, Size>::do_intersect(const Node &node, const Rectangle &rect) const
{
	if (node.is_leaf) {
		for (uint64_t iobj_ptr = node.leaf.obj_ptr_start_idx; iobj_ptr < node.leaf.obj_ptr_end_idx; ++iobj_ptr) {
			if (m_objs[m_obj_ptrs[iobj_ptr]].do_intersect(rect))
				return true;
		}
	} else {
		for (int iquad = 0; iquad < NUM_QUADS; ++iquad) {
			const Node &quad = m_nodes[node.node.kid_idx[iquad]];
			if (quad.stat.occupied_area > 0 && quad.arena.do_intersect(rect)) {
				if (quad.arena.is_inside(rect))
					return true;
				if (do_intersect(quad, rect))
					return true;
			}
		}
	}
	return false;
}

template <class T, class Size>
void StatQuadTree<T, Size>::intersect(const Rectangle &rect, Rectangles &intersection, vector<Size> &intersected_objs_indices) const
{
	if (m_intersected_objs.size() != m_objs.size())
		m_intersected_objs.resize(m_objs.size(), false);

	intersection.clear();
	intersected_objs_indices.clear();
	intersect(m_nodes.front(), rect, intersection, intersected_objs_indices);

	// m_intersected_objs was used to prevent multiple intersection with the same object; after intersection is done we need to reset it
	for (typename vector<Size>::const_iterator iobj = intersected_objs_indices.begin(); iobj != intersected_objs_indices.end(); ++iobj)
		m_intersected_objs[*iobj] = false;
}

template <class T, class Size>
void StatQuadTree<T, Size>::intersect(const Rectangle &rect, const DiagonalBand &band, Rectangles &intersection, vector<Size> &intersected_objs_indices) const
{
	if (band.do_contain(rect)) {
		intersect(rect, intersection, intersected_objs_indices);
		return;
	}

	if (m_intersected_objs.size() != m_objs.size())
		m_intersected_objs.resize(m_objs.size(), false);

	intersection.clear();
	intersected_objs_indices.clear();
	intersect(m_nodes.front(), rect, band, intersection, intersected_objs_indices);

	// m_intersected_objs was used to prevent multiple intersection with the same object; after intersection is done we need to reset it
	for (typename vector<Size>::const_iterator iobj = intersected_objs_indices.begin(); iobj != intersected_objs_indices.end(); ++iobj)
		m_intersected_objs[*iobj] = false;
}

template <class T, class Size>
void StatQuadTree<T, Size>::intersect(const Node &node, const Rectangle &rect, Rectangles &intersection, vector<Size> &intersected_objs_indices) const
{
	if (node.is_leaf) {
		for (uint64_t iobj_ptr = node.leaf.obj_ptr_start_idx; iobj_ptr < node.leaf.obj_ptr_end_idx; ++iobj_ptr) {
			Size obj_idx = m_obj_ptrs[iobj_ptr];
			if (!m_intersected_objs[obj_idx] && m_objs[obj_idx].do_intersect(rect)) {
				intersection.push_back(m_objs[obj_idx].intersect(rect));
				intersected_objs_indices.push_back(obj_idx);
				m_intersected_objs[obj_idx] = true;
			}
		}
	} else {
		for (int iquad = 0; iquad < NUM_QUADS; ++iquad) {
			const Node &quad = m_nodes[node.node.kid_idx[iquad]];
			if (quad.stat.occupied_area > 0 && quad.arena.do_intersect(rect))
				intersect(quad, rect, intersection, intersected_objs_indices);
		}
	}
}

template <class T, class Size>
void StatQuadTree<T, Size>::intersect(const Node &node, const Rectangle &rect, const DiagonalBand &band, Rectangles &intersection, vector<Size> &intersected_objs_indices) const
{
	if (node.is_leaf) {
		for (uint64_t iobj_ptr = node.leaf.obj_ptr_start_idx; iobj_ptr < node.leaf.obj_ptr_end_idx; ++iobj_ptr) {
			Size obj_idx = m_obj_ptrs[iobj_ptr];
			if (!m_intersected_objs[obj_idx] && m_objs[obj_idx].do_intersect(rect)) {
				Rectangle r(m_objs[obj_idx].intersect(rect));
				if (band.do_intersect(r)) {
					band.shrink2intersected(r);
					intersection.push_back(r);
					intersected_objs_indices.push_back(obj_idx);
					m_intersected_objs[obj_idx] = true;
				}
			}
		}
	} else {
		for (int iquad = 0; iquad < NUM_QUADS; ++iquad) {
			const Node &quad = m_nodes[node.node.kid_idx[iquad]];
			if (quad.stat.occupied_area > 0 && quad.arena.do_intersect(rect) && band.do_intersect(quad.arena.intersect(rect))) {
				// It is important to pass the original rect and not the shrinked one!!!
				// Otherwise the final objects will be intersected with the shrinked rect which was reduced by one of the quads.
				// The area of rect that intersects with other quads will be missed.
				// Also even if the the band contains the intersection between the quad and the rect we cannot omit the band.
				// Remember: at the end we can reach the object via only one the quads that contains it.
				intersect(quad, rect, band, intersection, intersected_objs_indices);
			}
		}
	}
}

template <class T, class Size>
void StatQuadTree<T, Size>::get_stat(const Rectangle &rect, Stat &stat) const
{
	stat.reset();
	get_stat(m_nodes.front(), rect, stat);
	if (!stat.occupied_area)
		stat.weighted_sum = stat.min_val = stat.max_val = numeric_limits<double>::quiet_NaN();
}

template <class T, class Size>
void StatQuadTree<T, Size>::get_stat(const Node &node, const Rectangle &rect, Stat &stat) const
{
	if (node.is_leaf) {
		for (uint64_t iobj_ptr = node.leaf.obj_ptr_start_idx; iobj_ptr < node.leaf.obj_ptr_end_idx; ++iobj_ptr) {
			const T &obj = m_objs[m_obj_ptrs[iobj_ptr]];
			Rectangle intersection = obj.intersect(rect, node.arena);
			if (intersection.is_non_empty_area())
				update_stat(obj, stat, intersection);
		}
	} else {
		for (int iquad = 0; iquad < NUM_QUADS; ++iquad) {
			const Node &quad = m_nodes[node.node.kid_idx[iquad]];
			if (quad.arena.do_intersect(rect)) {
				if (quad.arena.is_inside(rect)) {
					if (quad.stat.occupied_area)
						update_stat(quad, stat);
				} else
					get_stat(quad, rect, stat);
			}
		}
	}
}

template <class T, class Size>
void StatQuadTree<T, Size>::update_stat(const T &obj, Stat &stat, const Rectangle &intersection) const
{
	int64_t occupied_area = intersection.area();
	double val = obj.val(intersection, m_uptr);
	stat.weighted_sum += val * occupied_area;
	stat.min_val = min(val, stat.min_val);
	stat.max_val = max(val, stat.max_val);
	stat.occupied_area += occupied_area;
}

template <class T, class Size>
void StatQuadTree<T, Size>::update_stat(const Node &node, Stat &stat) const
{
	stat.weighted_sum += node.stat.weighted_sum;
	stat.min_val = min(node.stat.min_val, stat.min_val);
	stat.max_val = max(node.stat.max_val, stat.max_val);
	stat.occupied_area += node.stat.occupied_area;
}

template <class T, class Size>
void StatQuadTree<T, Size>::get_density_map(Rectangles &density_map) const
{
	density_map.clear();
	for (typename vector<Node>::const_iterator inode = m_nodes.begin(); inode != m_nodes.end(); ++inode) {
		if (inode->is_leaf)
			density_map.push_back(inode->arena);
	}
}

template <class T, class Size>
void StatQuadTree<T, Size>::serialize(BufferedFile &file)
{
	uint64_t size;

	file.write(&m_max_depth, sizeof(m_max_depth));
	file.write(&m_max_node_objs, sizeof(m_max_node_objs));

	size = m_nodes.size();
	file.write(&size, sizeof(size));

	size = m_obj_ptrs.size();
	file.write(&size, sizeof(size));

	size = m_obj_ptrs_free_chunks.size();
	file.write(&size, sizeof(size));

	size = m_objs.size();
	file.write(&size, sizeof(size));

	if (!m_nodes.empty())
		file.write(&m_nodes.front(), sizeof(m_nodes.front()) * m_nodes.size());

	if (!m_obj_ptrs.empty())
		file.write(&m_obj_ptrs.front(), sizeof(m_obj_ptrs.front()) * m_obj_ptrs.size());

	if (!m_obj_ptrs_free_chunks.empty())
		file.write(&m_obj_ptrs_free_chunks.front(), sizeof(m_obj_ptrs_free_chunks.front()) * m_obj_ptrs_free_chunks.size());

	if (!m_objs.empty())
		file.write(&m_objs.front(), sizeof(m_objs.front()) * m_objs.size());
}

template <class T, class Size>
void StatQuadTree<T, Size>::unserialize(BufferedFile &file)
{
	uint64_t size;

	file.read(&m_max_depth, sizeof(m_max_depth));
	file.read(&m_max_node_objs, sizeof(m_max_node_objs));

	file.read(&size, sizeof(size));
	m_nodes.resize(size);

	file.read(&size, sizeof(size));
	m_obj_ptrs.resize(size);

	file.read(&size, sizeof(size));
	m_obj_ptrs_free_chunks.resize(size);

	file.read(&size, sizeof(size));
	m_objs.resize(size);

	if (!m_nodes.empty())
		file.read(&m_nodes.front(), sizeof(m_nodes.front()) * m_nodes.size());

	if (!m_obj_ptrs.empty())
		file.read(&m_obj_ptrs.front(), sizeof(m_obj_ptrs.front()) * m_obj_ptrs.size());

	if (!m_obj_ptrs_free_chunks.empty())
		file.read(&m_obj_ptrs_free_chunks.front(), sizeof(m_obj_ptrs_free_chunks.front()) * m_obj_ptrs_free_chunks.size());

	if (!m_objs.empty())
		file.read(&m_objs.front(), sizeof(m_objs.front()) * m_objs.size());
}

template <class T, class Size>
void StatQuadTree<T, Size>::debug_print_tree() const
{
	debug_print_tree(m_nodes.front(), 0);
	printf("Objs: %ld\n", m_objs.size());
	printf("Nodes: %ld\n", m_nodes.size());
	printf("Obj ptrs: %ld\n", m_obj_ptrs.size());
	printf("Free chunks: %ld\n", m_obj_ptrs_free_chunks.size());
}

template <class T, class Size>
void StatQuadTree<T, Size>::debug_print_tree(const Node &node, unsigned depth) const
{
	printf("\n%*sArena: %s\n", depth * 2, "", node.arena.debug_str());
	printf("%*sIs leaf?: %d\n", (depth + 1) * 2, "", node.is_leaf);
	printf("%*sArea occupied: %ld\n", (depth + 1) * 2, "", node.stat.occupied_area);
	printf("%*sAvg: %g\tMin: %g\tMax: %g\n", (depth + 1) * 2, "", node.stat.occupied_area / (double)node.stat.weighted_sum, node.stat.min_val, node.stat.max_val);

	if (node.is_leaf) {
		printf("%*sKids: %ld, %ld, %ld\n", (depth + 1) * 2, "", node.leaf.obj_ptr_start_idx, node.leaf.obj_ptr_end_idx, node.leaf.obj_ptr_end_idx - node.leaf.obj_ptr_start_idx);
		for (uint64_t i = node.leaf.obj_ptr_start_idx; i < node.leaf.obj_ptr_end_idx; ++i) {
			printf("%*s%s", (depth + 2) * 2, "", m_objs[m_obj_ptrs[i]].debug_str());
			printf("\n");
		}
	} else {
		printf("%*sNW node\n", (depth + 1) * 2, "");
		debug_print_tree(m_nodes[node.node.kid_idx[NW]], depth + 1);
		printf("%*sNE node\n", (depth + 1) * 2, "");
		debug_print_tree(m_nodes[node.node.kid_idx[NE]], depth + 1);
		printf("%*sSE node\n", (depth + 1) * 2, "");
		debug_print_tree(m_nodes[node.node.kid_idx[SE]], depth + 1);
		printf("%*sSW node\n", (depth + 1) * 2, "");
		debug_print_tree(m_nodes[node.node.kid_idx[SW]], depth + 1);
	}
}

#endif /* QUADTREE_H_ */
