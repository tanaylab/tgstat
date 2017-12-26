#ifndef SEGMENT_FINDER_H_INCLUDED
#define SEGMENT_FINDER_H_INCLUDED

#include "Segment.h"

#include <queue>
#include <vector>

// SegmentsFinder stores a collection of segments and supports the following queries:
//     1. K nearest neighbors
//     2. Intersections and other queries, but not statistical queries (see: StatQuadTree) - TBD
//
// Note: the segments in SegmentsFinder can overlap.

template <class T>
class SegmentFinder {
public:
	typedef T ValueType;

	SegmentFinder() : m_root(NULL), m_max_depth(0), m_max_node_objs(0) { init(0, 0); }

	// start, end     - are the coordinates of the "arena" - the segment that should contain all the inserted objects.
	// max_depth      - maximal depth of the internal binary tree data structure.
	// max_node_objs  - maximal number of objects in the node. Once the number of objects exceed max_node_objs, the node is split.
	SegmentFinder(int64_t start, int64_t end, unsigned max_depth = 20, unsigned max_node_objs = 20);

	~SegmentFinder() { delete m_root; }

	void init(int64_t start, int64_t end, unsigned max_depth = 20, unsigned max_node_objs = 20);

	void reset();

	void reset(int64_t start, int64_t end);

	void insert(const T &obj);

	void debug_print_tree() const { debug_print_tree(m_root, 0); }

private:
#pragma pack(push)
#pragma pack(8)

	struct Node {
		Segment arena;
		vector<T> objs;
		Node *left;
		Node *right;

		Node() {}

		Node(const Segment &_arena) : arena(_arena), left(NULL), right(NULL) {}

		~Node() {
			delete left;
			delete right;
		}
	};

#pragma pack(pop)

	Node    *m_root;
	size_t   m_num_objs;
	unsigned m_max_depth;
	unsigned m_max_node_objs;

	void insert(Node *node, unsigned depth, const T &obj);

	void debug_print_tree(Node *node, unsigned depth) const;

public:
	// Nearest neighbor iterator returns the objects according to their proximity
	// to query segment starting from the nearest one.
	// Use excluded area to exclude objects that intersect that area.
	// Unlike standard STL iterators you can't copy this iterator nor compare it to another one.
	// Call is_end() to test whether the iterator reached the end.
	class NNIterator {
	public:
		NNIterator(SegmentFinder<T> *parent) : m_parent(parent) {}

		bool begin(const Segment &query, Segment excluded_area = Segment(0, -1)); // returns true if end is not reached yet
		bool next();                                                              // returns true if end is not reached yet
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
			bool operator<(const Neighbor &obj) const { return dist > obj.dist || dist == obj.dist && node; }
		};

		Segment                   m_query;
		Segment                   m_excluded_area;
		SegmentFinder<T>         *m_parent;
		priority_queue<Neighbor>  m_neighbors;

		T &cur_obj() { return *m_neighbors.top().obj; }
		const T &cur_obj() const { return *m_neighbors.top().obj; }

		void push_node(Node *node);

		// block copy constructor and default operators
		NNIterator(const NNIterator &) {}
		NNIterator &operator=(const NNIterator &) { return *this; }
		bool operator==(const NNIterator &) { return false; }
	};
};


//------------------------------ IMPLEMENTATION ----------------------------------------

//============================== NNIterator =============================================

template <class T>
bool SegmentFinder<T>::NNIterator::begin(const Segment &query, Segment excluded_area)
{
	m_query = query;
	m_excluded_area = excluded_area;
	m_neighbors = priority_queue<Neighbor>();

	if (!m_parent->m_num_objs)
		return false;

	push_node(m_parent->m_root);

	if (!m_neighbors.empty() && m_neighbors.top().obj)
		return true;

	return next();
}

template <class T>
bool SegmentFinder<T>::NNIterator::next()
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
		if (node->left) 
			push_node(node->left);
		if (node->right) 
			push_node(node->right);
	}
	return false;
}

template <class T>
void SegmentFinder<T>::NNIterator::push_node(Node *node)
{
	if (!m_excluded_area.do_contain(node->arena)) {
		m_neighbors.push(Neighbor(node, NULL, node->arena.dist2segment(m_query, true)));
		for (typename vector<T>::iterator iobj = node->objs.begin(); iobj != node->objs.end(); ++iobj) {
			if (!iobj->do_overlap(m_excluded_area))
				m_neighbors.push(Neighbor(NULL, &*iobj, iobj->dist2segment(m_query, true)));
		}
	}
}

//============================== SegmentFinder ==========================================

template <class T>
SegmentFinder<T>::SegmentFinder(int64_t start, int64_t end, unsigned max_depth, unsigned max_node_objs) :
	m_root(NULL),
	m_max_depth(0),
	m_max_node_objs(0)
{
	init(start, end, max_depth, max_node_objs);
}

template <class T>
void SegmentFinder<T>::init(int64_t start, int64_t end, unsigned max_depth, unsigned max_node_objs)
{
	m_max_depth = max_depth;
	m_max_node_objs = max_node_objs;

	reset(start, end);
}

template <class T>
void SegmentFinder<T>::reset()
{
	if (m_root) 
		reset(m_root->arena.start, m_root->arena.end);
}

template <class T>
void SegmentFinder<T>::reset(int64_t start, int64_t end)
{
	m_num_objs = 0;
	delete m_root;
	m_root = NULL;
	m_root = new Node(Segment(start, end));
}

template <class T>
void SegmentFinder<T>::insert(const T &obj)
{
	if (m_root->arena.do_overlap(obj))
		insert(m_root, 0, obj);
}

template <class T>
void SegmentFinder<T>::insert(Node *node, unsigned depth, const T &obj)
{
	// should the node be split?
	if (!node->left && node->objs.size() >= m_max_node_objs && depth < m_max_depth && node->arena.range() >= 2) {
		int64_t midpoint = (node->arena.start + node->arena.end) / 2;
		node->left = new Node(Segment(node->arena.start, midpoint));
		node->right = new Node(Segment(midpoint, node->arena.end));

		// reassign the objects to the newly created leaves
		typename vector<T>::iterator iobj = node->objs.begin();

		while (iobj < node->objs.end()) {
			if (node->left->arena.do_contain(*iobj)) {
				node->left->objs.push_back(*iobj);
				if (iobj < node->objs.end() - 1)
					*iobj = node->objs.back();
				node->objs.pop_back();    // vector::pop_back preserves the iterators
			} else if (node->right->arena.do_contain(*iobj)) {
				node->right->objs.push_back(*iobj);
				if (iobj < node->objs.end() - 1)
					*iobj = node->objs.back();
				node->objs.pop_back();    // vector::pop_back preserves the iterators
			} else
				++iobj;
		}
	}

	if (node->left && node->left->arena.do_contain(obj))
		insert(node->left, depth + 1, obj);
	else if (node->right && node->right->arena.do_contain(obj))
		insert(node->right, depth + 1, obj);
	else {
		node->objs.push_back(obj);
		m_num_objs++;
	}
}

template <class T>
void SegmentFinder<T>::debug_print_tree(Node *node, unsigned depth) const
{
	printf("\n%*sArena: %s\n", depth * 2, "", node->arena.debug_str());
	printf("%*sObjs: %ld\n", (depth + 1) * 2, "", node->objs.size());
	for (typename vector<T>::const_iterator iobj = node->objs.begin(); iobj < node->objs.end(); ++iobj) {
		printf("%*s%s", (depth + 2) * 2, "", iobj->debug_str());
		printf("\n");
	}

	if (node->left) {
		printf("%*sLEFT node\n", (depth + 1) * 2, "");
		debug_print_tree(node->left, depth + 1);
	}

	if (node->right) {
		printf("%*sRIGHT node\n", (depth + 1) * 2, "");
		debug_print_tree(node->right, depth + 1);
	}
}

#endif

