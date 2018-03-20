/*
 * StatQuadTreeCached.h
 *
 *  Created on: Jun 24, 2012
 *      Author: hoichman
 */

#ifndef STATQUADTREECACHED_H_
#define STATQUADTREECACHED_H_

#include <unordered_map>
#include <list>
#include <vector>

#include "DiagonalBand.h"
#include "StatQuadTree.h"
#include "TGLException.h"

// StatQuadTreeCached stores geographical objects and their value and allows efficient retrieval of various
// statistics about these objects.
//
// StatQuadTreeCached is constructed from an existing StatQuadTree. It is intended to supply efficient quering over huge number of objects
// that may or may not be fully loaded to the memory.
//
// StatQuadTreeCached is in many ways similar to StatQuadTree. However there are few important differences:
//
// 1. StatQuadTreeCached must be constructed and serialized prior to use of any "queries" such as get_stat, iterators, etc.
// 2. StatQuadTreeCached is constructed by calling serialize() function that accepts an instance of StatQuadTree and a file.
//    Indeed StatQuadTreeCached cannot exist without writing the whole tree to a file. StatQuadTreeCached also requires an existing StatQuadTree to
//    be constructed from.
// 3. Alternatively StatQuadTreeCached can be constructed and serialized via StatQuadTreeChachedSerializer - an adapter that allows also to
//    construct a tree with enormous number of objects. This can be handy when a StatQuadTree smply does not fit into memory.
// 4. StatQuadTreeCached does not support "insert" operations.
// 5. Unlike StatQuadTree that must be entirely stored in the memory, StatQuadTreeCached may load some parts of the quad tree
//    on demand in order to answer "get_stat" or "do_intersect" queries. Each part that StatQuadTreeCached loads is called a "chunk".
//    Chunks are stored in the memory and can be loaded or unloaded in accordance to their recent usage.
//    Chunk size and the maximal number of chunks StatQuadTreeCached stores at any given moment in the memory is controled by init() function.
// 6. Even on small trees StatQuadTreeCached may be more efficient than StatQuadTree since it stores objects in the memory in a more compact and
//    efficient way in the memo. This efficiency has a price: as mentioned earlier new objects cannot be inserted to StatQuadTreeCached.
// 7. StatQuadTreeCached replaces StatQuadTree::get_objs() function by an Iterator mechanism. The need to store in the memory all the objects
//    simultaneously is therefore avoided.

template <class T, class Size>
class StatQuadTreeCached {
public:
	template<class U, class V> friend class StatQuadTreeCachedSerializer;

	// nw = North West, ne = North East, se = South East, sw = South West;
	enum { NW, NE, SE, SW, NUM_QUADS };

	typedef T ValueType;

	typedef typename StatQuadTree<T, Size>::Stat Stat;

	// chunk size 0 means there's no caching at all
	StatQuadTreeCached() : m_chunk_size(0), m_max_num_chunks(0), m_num_chunks(0), m_bfile(NULL), m_uptr(NULL), m_local2global_id(NULL) {}
	StatQuadTreeCached(int64_t chunk_size, int64_t max_num_chunks) : m_uptr(NULL) { init(chunk_size, max_num_chunks); }

	~StatQuadTreeCached() { clear(); }

	// Subtrees whose size in bytes exceeds chunk_size will create a separate chunk.
	// max_num_chunks controls the maximal number of chunks stored in the memory at any given moment.
	// chunk_size == 0      => chunk size is unlimited
	// max_num_chunks == 0  => number of chunks in memory is unlimited
	void init(int64_t chunk_size, int64_t max_num_chunks);

	const Rectangle &get_arena() const { return m_root_chunk.top_node->arena; }

	// prior to serialize() file must be open for writing
	void serialize(BufferedFile &file, const StatQuadTree<T, Size> &qtree);

	// prior to unserialize() file must be open for reading;
	// BufferedFile object should remain undestructed until all queries are done
	void unserialize(BufferedFile &file);

	// resets all the data; unserialize must be called before next query
	void clear();

	// returns true if rect intersects any of the rectangles in the quad-tree
	bool do_intersect(const Rectangle &rect);

	// returns statitstics for the given rectangle: weighted sum, occupied area, minimal and maximal value of the object that intersect with rect.
	// Note: the average value can be achieved by dividing weighted sum by occupied.area
	void get_stat(const Rectangle &rect, Stat &stat);

	// Works in a similar way as a bandless get_stat. The band further restricts the query area.
	void get_stat(const Rectangle &rect, const DiagonalBand &band, Stat &stat);

	uint64_t get_num_objs() const { return m_num_objs; }

	bool empty() const { return !m_num_objs; }

	int64_t get_chunk_size() const { return m_chunk_size; }

	int64_t get_max_num_chunks() const { return m_max_num_chunks; }

	void set_uptr(void *uptr) { m_uptr = uptr; }

	void debug_print_tree();

private:
	// File format:
	// [num objs]
	// [root chunk start fpos]
	// <chunk 1>
	// <chunk 2>
	//   ...
	// <root chunk>
	//
	// Chunk format:
	//    [chunk size]
	//    [subtree top node offset from chunk start]
	//    NODES and LEAVES
	//
	// Node format:
	//    [struct Node]
	//
	// Leaf format:
	//    [struct Leaf]
	//    [struct Obj] x leaf.num_objs

#pragma pack(push)
#pragma pack(8)

	struct Obj {
		Size      id;
		T         obj;

		Obj(Size _id, const T &_obj) : id(_id), obj(_obj) {}
	};

	struct NodeBase {
		bool      is_leaf;
		Stat      stat;
		Rectangle arena;
	};

	struct Node : public NodeBase {
		int64_t kid_ptr[NUM_QUADS]; // positive values mark offset relative to the beginning of the chunk (nodes withing the chunk);
									// negative values mark the absolute file position of another chunk
	};

	struct Leaf : public NodeBase {
		unsigned  num_objs;
	};

#pragma pack(pop)

	struct Chunk {
		char     *mem;       // ptr to the memory block that holds the chunk
		NodeBase *top_node;  // ptr to the head node/leaf of the subtree
		int64_t   fpos;      // file position associated with the chunk

		Chunk() : mem(NULL), top_node(NULL) {}
	};

	class QuadRetriever {
	public:
		QuadRetriever(StatQuadTreeCached<T, Size> *parent, const Chunk &chunk, int64_t quad_ptr);
		~QuadRetriever();

		NodeBase    *quad() const { return m_quad; }
		const Chunk &chunk() const { return m_chunk; }
		int64_t      quad_ptr() const { return m_quad_ptr; }

	private:
		StatQuadTreeCached<T, Size> *m_parent;
		NodeBase *m_quad;
		Chunk     m_chunk;
		int64_t   m_quad_ptr;
	};

	typedef list<Chunk>                                       Chunks;
	typedef unordered_map<int64_t, typename Chunks::iterator> Fpos2ichunk;
	typedef vector<int64_t>                                   Stacked_chunks_fpos;

	int64_t             m_chunk_size;
	int64_t             m_max_num_chunks;
	int64_t             m_num_chunks;
	BufferedFile       *m_bfile;
	Chunk               m_root_chunk;
	Chunks              m_chunks;      // LRU list of chunks, the front element is the least recently used chunk
	Fpos2ichunk         m_fpos2ichunk;
	uint64_t            m_num_objs;
	vector<Size>       *m_local2global_id; // id map for tree objects (used for qtree serializationg as a chunk)

	// This preserves a chunk identified by fpos to be removed by LRU. Chunks that are in the stack are added here.
	// Since Iterator requires the stack to be preserved between the calls to next(), and other functions such as
	// get_stat() can be invoked meanwhile, it's not enough to save a list of chunks in the stack. We must keep a reference count too.
	Stacked_chunks_fpos m_stacked_chunks_fpos;
	void               *m_uptr;

	Obj *get_objs(NodeBase *node_base) { return (Obj *)((char *)node_base + sizeof(Leaf)); }

	// analyzes subtree and serializes it if needed; returns unserialized subtree size in bytes
	int64_t analyze_n_serialize_subtree(BufferedFile &file, const StatQuadTree<T, Size> &qtree, const typename StatQuadTree<T, Size>::Node &node, vector<int64_t> &chunks_fpos);

	// unconditionally serializes subtree, returns top node position relative to chunk_start_fpos
	int64_t serialize_subtree(BufferedFile &file, const StatQuadTree<T, Size> &qtree, const typename StatQuadTree<T, Size>::Node &node, vector<int64_t> &serialized_nodes_fpos, int64_t chunk_start_fpos);

	// serializes quad tree as a chunk (used to construct a huge tree from several subtree where each one is represented as a chunk), returns root chunk position
	int64_t serialize_as_chunk(BufferedFile &file, const StatQuadTree<T, Size> &qtree, vector<Size> &local2global_id);

	// gets a chunk residing at given file position
	const Chunk &get_chunk(int64_t fpos);
	void release_lru_chunk();
	Chunk &cache_chunk(const Chunk &chunk);

	bool do_intersect(const Chunk &chunk, const NodeBase *node_base, const Rectangle &rect);
	void get_stat(const Chunk &chunk, const NodeBase *node_base, const Rectangle &rect, Stat &stat);
	void get_stat(const Chunk &chunk, const NodeBase *node_base, const Rectangle &rect, const DiagonalBand &band, Stat &stat);
	void update_stat(const T &obj, Stat &stat, const Rectangle &intersection) const;
	void update_stat(const T &obj, Stat &stat, const Rectangle &intersection, const DiagonalBand &band) const;
	void update_stat(const NodeBase *node_base, Stat &stat) const;

	void debug_print_tree(const Chunk &chunk, NodeBase *node_base, unsigned depth);

public:
	// This iterator traces the tree. Unlike standard STL iterators you can't copy this iterator nor compare it to another one.
	// Call is_end() to test whether the iterator reached the end.
	class Iterator {
	public:
		Iterator(StatQuadTreeCached<T, Size> *parent) : m_parent(parent), m_cur_obj_idx(-1) {}
		~Iterator() { clear_stack(); }

		bool begin(); // returns true if end is not reached yet
		bool next();  // returns true if end is not reached yet
		bool is_end() const { return m_call_stack.empty(); }

		T &operator*() { return cur_obj(); }
		const T &operator*() const { return cur_obj(); }

		T *operator->() { return &cur_obj(); }
		const T *operator->() const { return &cur_obj(); }

		const Rectangle &containing_quad() const { return m_call_stack.back()->quad()->arena; }

	private:
		typedef vector<QuadRetriever *> CallStack;

		StatQuadTreeCached<T, Size> *m_parent;
		vector<bool>                 m_obj_used;    // STL implements it as a bit_vector
		CallStack                    m_call_stack;
		int                          m_cur_obj_idx;

		T &cur_obj() { return m_parent->get_objs(m_call_stack.back()->quad())[m_cur_obj_idx].obj; }
		const T &cur_obj() const { return m_parent->get_objs(m_call_stack.back()->quad())[m_cur_obj_idx].obj; }

		void clear_stack();

		// block copy constructor and default operators
		Iterator(const Iterator &) {}
		Iterator &operator=(const Iterator &) { return *this; }
		bool operator==(const Iterator &) { return false; }
	};

	// This iterator traces the tree and returns the leaves. Unlike standard STL iterators you can't copy this iterator nor compare it to another one.
	// Call is_end() to test whether the iterator reached the end.
	class DensityIterator {
	public:
		DensityIterator(StatQuadTreeCached<T, Size> *parent) : m_parent(parent) {}
		~DensityIterator() { clear_stack(); }

		bool begin(); // returns true if end is not reached yet
		bool next();  // returns true if end is not reached yet
		bool is_end() const { return m_call_stack.empty(); }

		Rectangle &operator*() { return cur_arena(); }
		const Rectangle &operator*() const { return cur_arena(); }

		Rectangle *operator->() { return &cur_arena(); }
		const Rectangle *operator->() const { return &cur_arena(); }

	private:
		typedef vector<QuadRetriever *> CallStack;

		StatQuadTreeCached<T, Size> *m_parent;
		CallStack                    m_call_stack;

		Rectangle &cur_arena() { return m_call_stack.back()->quad()->arena; }
		const Rectangle &cur_arena() const { return m_call_stack.back()->quad()->arena; }

		void clear_stack();

		// block copy constructor and default operators
		DensityIterator(const DensityIterator &) {}
		DensityIterator &operator=(const DensityIterator &) { return *this; }
		bool operator==(const DensityIterator &) { return false; }
	};
};

typedef StatQuadTreeCached<Rectangle_val<float>, uint64_t> RectsQuadTreeCached;
typedef StatQuadTreeCached<Point_val<float>, uint64_t>     PointsQuadTreeCached;

//------------------------------ IMPLEMENTATION ----------------------------------------

//================= QuadRetriever ==================

template <class T, class Size>
StatQuadTreeCached<T, Size>::QuadRetriever::QuadRetriever(StatQuadTreeCached<T, Size> *parent, const Chunk &chunk, int64_t quad_ptr) :
	m_parent(parent)
{
	m_quad_ptr = quad_ptr;

	if (quad_ptr > 0) {    // quad_ptr is a memory offset relative to the beginning of chunk
		m_quad = (NodeBase *)(chunk.mem + m_quad_ptr);
		m_chunk = chunk;
	} else {               // quad_ptr is a file position of the chunk
		m_chunk = m_parent->get_chunk(m_quad_ptr);
		m_quad = m_chunk.top_node;
		m_parent->m_stacked_chunks_fpos.push_back(m_quad_ptr);
	}
}

template <class T, class Size>
StatQuadTreeCached<T, Size>::QuadRetriever::~QuadRetriever()
{
	if (m_quad_ptr < 0) {
		// most of the chances that the chunk is at the end of the list, but it might not be the case if two Iterators are used
		if (m_parent->m_stacked_chunks_fpos.back() != m_quad_ptr) {
			Stacked_chunks_fpos::iterator itr = find(m_parent->m_stacked_chunks_fpos.begin(), m_parent->m_stacked_chunks_fpos.end(), m_quad_ptr);
			*itr = m_parent->m_stacked_chunks_fpos.back();
		}
		m_parent->m_stacked_chunks_fpos.pop_back();
	}
}

//================= Iterator ==================

template <class T, class Size>
void StatQuadTreeCached<T, Size>::Iterator::clear_stack()
{
	for (typename CallStack::reverse_iterator icall = m_call_stack.rbegin(); icall != m_call_stack.rend(); ++icall)
		delete *icall;
}

template <class T, class Size>
bool StatQuadTreeCached<T, Size>::Iterator::begin()
{
	clear_stack();
	m_call_stack.push_back(new QuadRetriever(m_parent, m_parent->m_root_chunk, m_parent->m_root_chunk.fpos));
	m_cur_obj_idx = -1;
	m_obj_used.resize(m_parent->get_num_objs(), 0);
	return next();
}

template <class T, class Size>
bool StatQuadTreeCached<T, Size>::Iterator::next()
{
	int64_t last_quad_ptr = 0;

	++m_cur_obj_idx;

	while (1) {
		if (m_call_stack.empty())
			return false;

		QuadRetriever *qr = m_call_stack.back();

		if (qr->quad()->is_leaf) {
			Leaf *leaf = (Leaf *)qr->quad();

			if ((unsigned)m_cur_obj_idx < leaf->num_objs) {
				Size id = m_parent->get_objs(m_call_stack.back()->quad())[m_cur_obj_idx].id;

				if (m_obj_used[id]) {
					++m_cur_obj_idx;
					continue;
				}

				m_obj_used[id] = true;
				return true;
			}

			m_cur_obj_idx = 0;
			last_quad_ptr = qr->quad_ptr();
			delete qr;
			m_call_stack.pop_back();
			continue;
		}

		Node *node = (Node *)qr->quad();

		// tracing up the tree
		if (last_quad_ptr) {
			// all the quads had been explored => go back
			if (last_quad_ptr == node->kid_ptr[NUM_QUADS - 1]) {
				last_quad_ptr = qr->quad_ptr();
				delete qr;
				m_call_stack.pop_back();
				continue;
			}

			for (int iquad = 0; iquad < NUM_QUADS - 1; ++iquad) {
				if (last_quad_ptr == node->kid_ptr[iquad]) {
					// trace down the next kid
					m_call_stack.push_back(new QuadRetriever(m_parent, qr->chunk(), node->kid_ptr[iquad + 1]));
					last_quad_ptr = 0;
					break;
				}
			}
		} else      // trace down the first kid
			m_call_stack.push_back(new QuadRetriever(m_parent, qr->chunk(), node->kid_ptr[0]));
	}
}

//================= DensityIterator ==================

template <class T, class Size>
void StatQuadTreeCached<T, Size>::DensityIterator::clear_stack()
{
	for (typename CallStack::reverse_iterator icall = m_call_stack.rbegin(); icall != m_call_stack.rend(); ++icall)
		delete *icall;
}

template <class T, class Size>
bool StatQuadTreeCached<T, Size>::DensityIterator::begin()
{
	clear_stack();
	m_call_stack.push_back(new QuadRetriever(m_parent, m_parent->m_root_chunk, m_parent->m_root_chunk.fpos));

	if (m_call_stack.back()->quad()->is_leaf)
		return true;

	return next();
}

template <class T, class Size>
bool StatQuadTreeCached<T, Size>::DensityIterator::next()
{
	int64_t last_quad_ptr = 0;

	while (1) {
		if (m_call_stack.empty())
			return false;

		QuadRetriever *qr = m_call_stack.back();

		if (qr->quad()->is_leaf) {
			last_quad_ptr = qr->quad_ptr();
			delete qr;
			m_call_stack.pop_back();
			continue;
		}

		Node *node = (Node *)qr->quad();

		// tracing up the tree
		if (last_quad_ptr) {
			// all the quads had been explored => go back
			if (last_quad_ptr == node->kid_ptr[NUM_QUADS - 1]) {
				last_quad_ptr = qr->quad_ptr();
				delete qr;
				m_call_stack.pop_back();
				continue;
			}

			for (int iquad = 0; iquad < NUM_QUADS - 1; ++iquad) {
				if (last_quad_ptr == node->kid_ptr[iquad]) {
					// trace down the next kid
					m_call_stack.push_back(new QuadRetriever(m_parent, qr->chunk(), node->kid_ptr[iquad + 1]));
					if (m_call_stack.back()->quad()->is_leaf)
						return true;

					last_quad_ptr = 0;
					break;
				}
			}
		} else {     // trace down the first kid
			m_call_stack.push_back(new QuadRetriever(m_parent, qr->chunk(), node->kid_ptr[0]));
			if (m_call_stack.back()->quad()->is_leaf)
				return true;
		}
	}
}

//================= StatQuadTreeCached ==================

template <class T, class Size>
void StatQuadTreeCached<T, Size>::init(int64_t chunk_size, int64_t max_num_chunks)
{
	m_chunk_size = chunk_size;
	m_max_num_chunks = max_num_chunks;
	clear();
}

template <class T, class Size>
void StatQuadTreeCached<T, Size>::serialize(BufferedFile &file, const StatQuadTree<T, Size> &qtree)
{
	m_num_objs = qtree.get_num_objs();
	file.write(&m_num_objs, sizeof(m_num_objs));

	if (m_num_objs) {
		int64_t tree_start_fpos = file.tell();
		int64_t root_chunk_start_fpos = 0;
		vector<int64_t> serialized_nodes_fpos(qtree.m_nodes.size(), 0);

		file.write(&root_chunk_start_fpos, sizeof(root_chunk_start_fpos)); // write a placeholder for the file position of the root chunk
		analyze_n_serialize_subtree(file, qtree, qtree.m_nodes.front(), serialized_nodes_fpos);

		int64_t cur_fpos = file.tell();
		file.seek(tree_start_fpos, SEEK_SET);
		root_chunk_start_fpos = serialized_nodes_fpos.front();
		file.write(&root_chunk_start_fpos, sizeof(root_chunk_start_fpos));
		file.seek(cur_fpos, SEEK_SET);
	}

	if (file.error())
		TGLError< StatQuadTreeCached<T, Size> >("Writing file %s: %s", file.file_name().c_str(), strerror(errno));
}

template <class T, class Size>
int64_t StatQuadTreeCached<T, Size>::serialize_as_chunk(BufferedFile &file, const StatQuadTree<T, Size> &qtree, vector<Size> &local2global_id)
{
	vector<int64_t> serialized_nodes_fpos(qtree.m_nodes.size(), 0);

	m_local2global_id = &local2global_id;
	analyze_n_serialize_subtree(file, qtree, qtree.m_nodes.front(), serialized_nodes_fpos);
	m_local2global_id = NULL;
	return serialized_nodes_fpos.front();
}

template <class T, class Size>
int64_t StatQuadTreeCached<T, Size>::analyze_n_serialize_subtree(BufferedFile &file, const StatQuadTree<T, Size> &qtree, const typename StatQuadTree<T, Size>::Node &node, vector<int64_t> &chunks_fpos)
{
	int64_t size = 0;

	if (node.is_leaf)
		size = sizeof(Leaf) + (node.leaf.obj_ptr_end_idx - node.leaf.obj_ptr_start_idx) * sizeof(Obj);
	else {
		for (int iquad = 0; iquad < NUM_QUADS; iquad++) {
			unsigned quad_idx = node.node.kid_idx[iquad];
			int64_t subtree_size = analyze_n_serialize_subtree(file, qtree, qtree.m_nodes[quad_idx], chunks_fpos);

			if (subtree_size)
				size += subtree_size;
		}
		size += sizeof(Node);
	}

	if ((m_chunk_size && size > m_chunk_size) || &node == &qtree.m_nodes.front()) {
		int64_t chunk_start_fpos = file.tell();
		int64_t top_node_offset = 0;

		size += sizeof(size) + sizeof(top_node_offset);
		file.write(&size, sizeof(size));                       // write the chunk size
		file.write(&top_node_offset, sizeof(top_node_offset)); // write a placeholder for a top node file position

		top_node_offset = serialize_subtree(file, qtree, node, chunks_fpos, chunk_start_fpos);

		int64_t cur_fpos = file.tell();
		file.seek(chunk_start_fpos + sizeof(size), SEEK_SET);
		file.write(&top_node_offset, sizeof(top_node_offset));
		file.seek(cur_fpos, SEEK_SET);

		chunks_fpos[&node - &qtree.m_nodes.front()] = chunk_start_fpos;
		return 0;
	}

	return size;
}

template <class T, class Size>
int64_t StatQuadTreeCached<T, Size>::serialize_subtree(BufferedFile &file, const StatQuadTree<T, Size> &qtree, const typename StatQuadTree<T, Size>::Node &node, vector<int64_t> &chunks_fpos, int64_t chunk_start_fpos)
{
	int64_t top_node_fpos;

	if (node.is_leaf) {
		Leaf new_leaf;

		// ensure that stucture's members are padded with zeroes (and not some junk) otherwise there will be no binary consistency between the files with identical quad trees
		memset(&new_leaf, 0, sizeof(new_leaf));
		new_leaf.is_leaf = true;
		new_leaf.stat = node.stat;
		new_leaf.arena = node.arena;
		new_leaf.num_objs = node.leaf.obj_ptr_end_idx - node.leaf.obj_ptr_start_idx;
		top_node_fpos = file.tell();
		file.write(&new_leaf, sizeof(new_leaf));

		for (uint64_t iobj_ptr = node.leaf.obj_ptr_start_idx; iobj_ptr < node.leaf.obj_ptr_end_idx; ++iobj_ptr) {
			Size obj_ptr = qtree.m_obj_ptrs[iobj_ptr];
			Size id = m_local2global_id ? (*m_local2global_id)[obj_ptr] : obj_ptr;

			// Obj is a structure. It's members might be padded.
			// Therefore writing each of the structure members separately might not give the same result as writing the whole structure as a bulk.
			// If this happens get_objs() will incorrectly make the casting to Obj *, and which might result in memory corruption.
			// Therefore we create an instance of Obj, fill its members and write the structure as a bulk.
			Obj obj(id, qtree.m_objs[obj_ptr]);
			file.write(&obj, sizeof(obj));
		}
	} else {
		Node new_node;

		// ensure that stucture's members are padded with zeroes (and not some junk) otherwise there will be no binary consistency between the files with identical quad trees
		memset(&new_node, 0, sizeof(new_node));
		new_node.is_leaf = false;
		new_node.stat = node.stat;
		new_node.arena = node.arena;

		for (int iquad = 0; iquad < NUM_QUADS; iquad++) {
			unsigned quad_idx = node.node.kid_idx[iquad];

			new_node.kid_ptr[iquad] = -chunks_fpos[quad_idx]; // file position is negative to distinguish it from memory offset

			if (!new_node.kid_ptr[iquad]) // quad does not point to another chunk => write the quad
				// that will be the memory offset relative to the beginning of chunk
				new_node.kid_ptr[iquad] = serialize_subtree(file, qtree, qtree.m_nodes[quad_idx], chunks_fpos, chunk_start_fpos);
		}

		top_node_fpos = file.tell();
		file.write(&new_node, sizeof(new_node));
	}
	return top_node_fpos - chunk_start_fpos;
}

template <class T, class Size>
void StatQuadTreeCached<T, Size>::unserialize(BufferedFile &file)
{
	clear();
	m_bfile = &file;

	if (file.read(&m_num_objs, sizeof(m_num_objs)) != sizeof(m_num_objs)) {
		if (m_bfile->error())
			TGLError< StatQuadTreeCached<T, Size> >("Reading file %s: %s", m_bfile->file_name().c_str(), strerror(errno));
		TGLError< StatQuadTreeCached<T, Size> >("Invalid format of file %s", m_bfile->file_name().c_str());
	}

	if (m_num_objs) {
		int64_t root_chunk_start_fpos;

		if (file.read(&root_chunk_start_fpos, sizeof(root_chunk_start_fpos)) != sizeof(root_chunk_start_fpos)) {
			if (m_bfile->error())
				TGLError< StatQuadTreeCached<T, Size> >("Reading file %s: %s", m_bfile->file_name().c_str(), strerror(errno));
			TGLError< StatQuadTreeCached<T, Size> >("Invalid format of file %s", m_bfile->file_name().c_str());
		}

		m_root_chunk = get_chunk(-root_chunk_start_fpos);
		m_stacked_chunks_fpos.push_back(m_root_chunk.fpos); // never let the root chunk to be unloaded from LRU list
	}
}

template <class T, class Size>
const typename StatQuadTreeCached<T, Size>::Chunk &StatQuadTreeCached<T, Size>::get_chunk(int64_t fpos)
{
	typename Fpos2ichunk::iterator ifpos2ichunk = m_fpos2ichunk.find(fpos);

	if (ifpos2ichunk != m_fpos2ichunk.end()) {
		typename Chunks::iterator ichunk = ifpos2ichunk->second;
		// move the chunk to the end of LRU list
		m_chunks.push_front(*ichunk);
		m_chunks.erase(ichunk);
		ichunk = m_chunks.begin();
		ifpos2ichunk->second = ichunk;
		return *ichunk;
	}

	release_lru_chunk();

	int64_t size;
	Chunk chunk;

	m_bfile->seek(-fpos, SEEK_SET);   // chunk fpos is always a negative number to distinguish it from memory offset

	if (m_bfile->read(&size, sizeof(size)) != sizeof(size)) {
		if (m_bfile->error())
			TGLError< StatQuadTreeCached<T, Size> >("Reading file %s: %s", m_bfile->file_name().c_str(), strerror(errno));
		TGLError< StatQuadTreeCached<T, Size> >("Invalid format of file %s", m_bfile->file_name().c_str());
	}

	chunk.mem = new char[size];
	*(int64_t *)chunk.mem = size;
	if (m_bfile->read(chunk.mem + sizeof(size), size - sizeof(size)) != size - sizeof(size)) {
		if (m_bfile->error())
			TGLError< StatQuadTreeCached<T, Size> >("Reading file %s: %s", m_bfile->file_name().c_str(), strerror(errno));
		TGLError< StatQuadTreeCached<T, Size> >("Invalid format of file %s", m_bfile->file_name().c_str());
	}

	int64_t top_node_offset = *(int64_t *)(chunk.mem + sizeof(size));
	chunk.top_node = (NodeBase *)(chunk.mem + top_node_offset);
	chunk.fpos = fpos;

	return cache_chunk(chunk);
}

template <class T, class Size>
void StatQuadTreeCached<T, Size>::release_lru_chunk()
{
	// release least recenetly used chunk if needed
	if (m_max_num_chunks && m_num_chunks >= m_max_num_chunks) {
		for (typename Chunks::reverse_iterator ichunk = m_chunks.rbegin(); ichunk != m_chunks.rend(); ++ichunk) {
			// never remove chunks that are currently in the calling stack
			if (find(m_stacked_chunks_fpos.begin(), m_stacked_chunks_fpos.end(), ichunk->fpos) == m_stacked_chunks_fpos.end()) {
				m_fpos2ichunk.erase(ichunk->fpos);
				delete []ichunk->mem;
				m_chunks.erase((++ichunk).base());
				--m_num_chunks;
				return;
			}
		}
	}
}

template <class T, class Size>
typename StatQuadTreeCached<T, Size>::Chunk &StatQuadTreeCached<T, Size>::cache_chunk(const StatQuadTreeCached<T, Size>::Chunk &chunk)
{
	m_chunks.push_front(chunk);
	++m_num_chunks;
	m_fpos2ichunk[chunk.fpos] = m_chunks.begin();
	return m_chunks.front();
}

template <class T, class Size>
void StatQuadTreeCached<T, Size>::clear()
{
	// don't delete m_bfile, we're not the owner of the pointer
	m_bfile = NULL;

	m_fpos2ichunk.clear();
	for (typename Chunks::iterator ichunk = m_chunks.begin(); ichunk != m_chunks.end(); ++ichunk)
		delete []ichunk->mem;
	m_chunks.clear();
	m_num_chunks = 0;
	m_stacked_chunks_fpos.clear();
	m_root_chunk = Chunk();
	m_num_objs = 0;
	m_local2global_id = NULL;
}


template <class T, class Size>
bool StatQuadTreeCached<T, Size>::do_intersect(const Rectangle &rect)
{
	if (empty())
		return false;

	return do_intersect(m_root_chunk, m_root_chunk.top_node, rect);
}

template <class T, class Size>
bool StatQuadTreeCached<T, Size>::do_intersect(const Chunk &chunk, const NodeBase *node_base, const Rectangle &rect)
{
	if (node_base->is_leaf) {
		Leaf *leaf = (Leaf *)node_base;
		// objects start right after num_objs member (yeah, the arithmetics is quite ugly, I know...)
		Obj *objs = get_objs(leaf);

		for (unsigned i = 0; i < leaf->num_objs; ++i) {
			const T &obj = objs[i].obj;
			if (obj.do_intersect(rect))
				return true;
		}
	} else {
		Node *node = (Node *)node_base;
		for (int iquad = 0; iquad < NUM_QUADS; ++iquad) {
			QuadRetriever qr(this, chunk, node->kid_ptr[iquad]);
			if (qr.quad()->stat.occupied_area > 0 && qr.quad()->arena.do_intersect(rect)) {
				if (qr.quad()->arena.is_inside(rect))
					return true;
				if (do_intersect(qr.chunk(), qr.quad(), rect))
					return true;
			}
		}
	}
	return false;
}

template <class T, class Size>
void StatQuadTreeCached<T, Size>::get_stat(const Rectangle &rect, Stat &stat)
{
	stat.reset();

	if (!empty())
		get_stat(m_root_chunk, m_root_chunk.top_node, rect, stat);

	if (!stat.occupied_area)
		stat.weighted_sum = stat.min_val = stat.max_val = numeric_limits<double>::quiet_NaN();
}

template <class T, class Size>
void StatQuadTreeCached<T, Size>::get_stat(const Rectangle &rect, const DiagonalBand &band, Stat &stat)
{
	stat.reset();

	if (!empty())
		get_stat(m_root_chunk, m_root_chunk.top_node, rect, band, stat);

	if (!stat.occupied_area)
		stat.weighted_sum = stat.min_val = stat.max_val = numeric_limits<double>::quiet_NaN();
}

template <class T, class Size>
void StatQuadTreeCached<T, Size>::get_stat(const Chunk &chunk, const NodeBase *node_base, const Rectangle &rect, Stat &stat)
{
	if (node_base->is_leaf) {
		Leaf *leaf = (Leaf *)node_base;
		// objects start right after num_objs member (yeah, the arithmetics is quite ugly, I know...)
		Obj *objs = get_objs(leaf);

		for (unsigned i = 0; i < leaf->num_objs; ++i) {
			const T &obj = objs[i].obj;
			Rectangle intersection = obj.intersect(rect, leaf->arena);
			if (intersection.is_non_empty_area())
				update_stat(obj, stat, intersection);
		}
	} else {
		Node *node = (Node *)node_base;
		for (int iquad = 0; iquad < NUM_QUADS; ++iquad) {
			QuadRetriever qr(this, chunk, node->kid_ptr[iquad]);
			if (qr.quad()->arena.do_intersect(rect)) {
				if (qr.quad()->arena.is_inside(rect)) {
					if (qr.quad()->stat.occupied_area)
						update_stat(qr.quad(), stat);
				} else
					get_stat(qr.chunk(), qr.quad(), rect, stat);
			}
		}
	}
}

template <class T, class Size>
void StatQuadTreeCached<T, Size>::get_stat(const Chunk &chunk, const NodeBase *node_base, const Rectangle &rect, const DiagonalBand &band, Stat &stat)
{
	if (node_base->is_leaf) {
		Leaf *leaf = (Leaf *)node_base;
		// objects start right after num_objs member (yeah, the arithmetics is quite ugly, I know...)
		Obj *objs = get_objs(leaf);

		for (unsigned i = 0; i < leaf->num_objs; ++i) {
			const T &obj = objs[i].obj;
			Rectangle intersection = obj.intersect(rect, leaf->arena);
			if (intersection.is_non_empty_area()) {
				if (band.do_contain(intersection))
					update_stat(obj, stat, intersection);
				else if (band.do_intersect(intersection)) {
					band.shrink2intersected(intersection);
					update_stat(obj, stat, intersection, band);
				}
			}
		}
	} else {
		Node *node = (Node *)node_base;
		for (int iquad = 0; iquad < NUM_QUADS; ++iquad) {
			QuadRetriever qr(this, chunk, node->kid_ptr[iquad]);
			if (qr.quad()->arena.do_intersect(rect)) {
				if (qr.quad()->arena.is_inside(rect)) {
					if (band.do_contain(qr.quad()->arena))
						update_stat(qr.quad(), stat);
					else if (band.do_intersect(qr.quad()->arena)) {
						Rectangle r(qr.quad()->arena);
						band.shrink2intersected(r);
						get_stat(qr.chunk(), qr.quad(), r, band, stat);
					}
				} else {
					Rectangle intersection = qr.quad()->arena.intersect(rect);
					if (band.do_contain(intersection))
						get_stat(qr.chunk(), qr.quad(), rect, stat);
					else if (band.do_intersect(intersection)) {
						Rectangle r(qr.quad()->arena);
						band.shrink2intersected(r);
						get_stat(qr.chunk(), qr.quad(), intersection, band, stat);
					}
				}
			}
		}
	}
}

template <class T, class Size>
void StatQuadTreeCached<T, Size>::update_stat(const T &obj, Stat &stat, const Rectangle &intersection) const
{
	int64_t occupied_area = intersection.area();
	double val = obj.val(intersection, m_uptr);
	stat.weighted_sum += val * occupied_area;
	stat.min_val = min(val, stat.min_val);
	stat.max_val = max(val, stat.max_val);
	stat.occupied_area += occupied_area;
}


template <class T, class Size>
void StatQuadTreeCached<T, Size>::update_stat(const T &obj, Stat &stat, const Rectangle &intersection, const DiagonalBand &band) const
{
	int64_t occupied_area = band.intersected_area(intersection);
	double val = obj.val(intersection, band, m_uptr);
	stat.weighted_sum += val * occupied_area;
	stat.min_val = min(val, stat.min_val);
	stat.max_val = max(val, stat.max_val);
	stat.occupied_area += occupied_area;
}

template <class T, class Size>
void StatQuadTreeCached<T, Size>::update_stat(const NodeBase *node_base, Stat &stat) const
{
	stat.weighted_sum += node_base->stat.weighted_sum;
	stat.min_val = min(node_base->stat.min_val, stat.min_val);
	stat.max_val = max(node_base->stat.max_val, stat.max_val);
	stat.occupied_area += node_base->stat.occupied_area;
}

template <class T, class Size>
void StatQuadTreeCached<T, Size>::debug_print_tree()
{
	if (!empty())
		debug_print_tree(m_root_chunk, m_root_chunk.top_node, 0);
	printf("Objs: %llu\n", m_num_objs);
}

template <class T, class Size>
void StatQuadTreeCached<T, Size>::debug_print_tree(const Chunk &chunk, NodeBase *node_base, unsigned depth)
{
	Rectangle arena;
	printf("\n%*sArena: %s\n", depth * 2, "", node_base->arena.debug_str());
	printf("%*sIs leaf?: %d\n", (depth + 1) * 2, "", node_base->is_leaf);
	printf("%*sArea occupied: %lld\n", (depth + 1) * 2, "", node_base->stat.occupied_area);
	printf("%*sAvg: %g\tMin: %g\tMax: %g\n", (depth + 1) * 2, "", node_base->stat.occupied_area / (double)node_base->stat.weighted_sum, node_base->stat.min_val, node_base->stat.max_val);

	if (node_base->is_leaf) {
		Leaf *leaf = (Leaf *)node_base;
		Obj *objs = get_objs(leaf);

		printf("%*sKids: %d\n", (depth + 1) * 2, "", leaf->num_objs);
		for (unsigned i = 0; i < leaf->num_objs; ++i) {
			printf("%*s%s", (depth + 2) * 2, "", objs[i].obj.debug_str());
			printf("\n");
		}
	} else {
		Node *node = (Node *)node_base;
		for (int iquad = 0; iquad < NUM_QUADS; ++iquad) {
			QuadRetriever qr(this, chunk, node->kid_ptr[iquad]);
			if (iquad == NW)
				printf("%*sNW node\n", (depth + 1) * 2, "");
			else if (iquad == NE)
				printf("%*sNE node\n", (depth + 1) * 2, "");
			else if (iquad == SW)
				printf("%*sSW node\n", (depth + 1) * 2, "");
			else
				printf("%*sSE node\n", (depth + 1) * 2, "");
			debug_print_tree(qr.chunk(), qr.quad(), depth + 1);
		}
	}
}

#endif /* STATQUADTREECACHED_H_ */
