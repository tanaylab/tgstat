#ifndef stdalg_ds_Matrix_h
#define stdalg_ds_Matrix_h 1

#include <vector>

template<class T>
class MatrixColIter;
template<class T>
bool operator!=(const MatrixColIter<T> &i1, const MatrixColIter<T> &i2);
template<class T>
class MatrixColIter {

public:

	friend bool operator!=<T>(const MatrixColIter<T> &i1, const MatrixColIter<T> &i2);

private:

	typedef typename vector<T>::const_iterator vec_iter;
	vec_iter m_pos;

	int m_row_size;

public:


	MatrixColIter(vec_iter pos, int row_size) :
		m_pos(pos),
		m_row_size(row_size)
	{}

	void operator++() {
		m_pos += m_row_size;
	}
	void operator--() {
		m_pos -= m_row_size;
	}

	T operator*() {
		return(*m_pos);
	}
};

template<class T>
bool operator!=(const MatrixColIter<T> &i1, const MatrixColIter<T> &i2) {
	return(i1.m_pos != i2.m_pos);
}

template<class T>
class ConstMatrixRow {

private:
	typedef typename vector<T>::const_iterator vec_citer;
	vec_citer m_pos;

public:

	ConstMatrixRow(vec_citer pos) :
		m_pos(pos)
	{}

	const T &operator[](int col) const { return(*(m_pos + col)); }
};

template<class T>
class MatrixRow {

private:
	typedef typename vector<T>::iterator vec_iter;
	vec_iter m_pos;

public:

	MatrixRow(vec_iter pos) :
		m_pos(pos)
	{}

	T &operator[](int col) {
		return(*(m_pos + col));
	}
	const T &operator[](int col) const { return(*(m_pos + col)); }
};

template<class T>
class ConstMatrixCol {

private:
	friend class MatrixColIter<T>;

	typedef typename vector<T>::const_iterator vec_iter;
	vec_iter m_pos;

	int m_row_size;
	int m_col_size;

public:

	ConstMatrixCol(vec_iter pos, int row_size, int col_size) :
		m_pos(pos),
		m_row_size(row_size),
		m_col_size(col_size)
	{}

	const T &operator[](int row) const {
		return(*(m_pos + row * m_row_size));
       	}
	MatrixColIter<T> begin() {
		return(MatrixColIter<T>(m_pos, m_row_size));
	}

	MatrixColIter<T> end() {
		return(MatrixColIter<T>(m_pos + m_row_size * m_col_size, m_row_size));
	}

};

template<class T>
class MatrixCol {

private:
	friend class MatrixColIter<T>;

	typedef typename vector<T>::iterator vec_iter;
	vec_iter m_pos;

	int m_row_size;
	int m_col_size;

public:

	MatrixCol(vec_iter pos, int row_size, int col_size) :
		m_pos(pos),
		m_row_size(row_size),
		m_col_size(col_size)
	{}

	T &operator[](int row) {
		return(*(m_pos + row * m_row_size));
	}
	const T &operator[](int row) const { return(*(m_pos + row * m_row_size)); }
	MatrixColIter<T> begin() {
		return(MatrixColIter<T>(m_pos, m_row_size));
	}

	MatrixColIter<T> end() {
		return(MatrixColIter<T>(m_pos + m_row_size * m_col_size, m_row_size));
	}

};


template<class T>
class Matrix {
protected:

	typedef typename vector<T>::iterator vec_iter;
	vector<T> m_mat;

	int m_num_cols;
	int m_num_rows;

public:

	int col_size() const {
		return(m_num_cols);
	}
	int row_size() const {
		return(m_num_rows);
	}

	Matrix() :
		m_num_cols(0),
		m_num_rows(0)
	{}

	Matrix(int num_rows, int num_cols, T def = T()) :
		m_mat(num_cols * num_rows, def),
		m_num_cols(num_cols),
		m_num_rows(num_rows)
	{
//		resize(num_rows, num_cols, def);
	}

	void resize_rows(int num_rows, T def = T()) {
		m_num_rows = num_rows;
		m_mat.resize(num_rows * m_num_cols, def);
	}
	void resize(int num_rows, int num_cols, T def = T());

	vector<T> &get_vector() {
		return(m_mat);
	}
	const vector<T> &get_vector() const {
		return(m_mat);
	}

	T get_elem(int row_id, int col_id) const {
		return(m_mat[row_id * m_num_cols + col_id]);
	}
	MatrixRow<T> operator[](int row_id) {
		MatrixRow<T> row(m_mat.begin() + row_id * m_num_cols);
		return(row);
	}
	ConstMatrixRow<T> operator[](int row_id) const {
		ConstMatrixRow<T> row(m_mat.begin() + row_id * m_num_cols);
		return(row);
	}
	MatrixCol<T> get_col(int col_id) {
		MatrixCol<T> col(m_mat.begin() + col_id, m_num_cols, m_num_rows);
		return(col);
	}
	ConstMatrixCol<T> get_col(int col_id) const {
		ConstMatrixCol<T> col(m_mat.begin() + col_id, m_num_cols, m_num_rows);
		return(col);
	}

	void fill(T dup) {
		std::fill(m_mat.begin(), m_mat.end(), dup);
	}
	
	void get_column_vector(int col, vector<T> &vec) const {
		vec.resize(0);
		ConstMatrixCol<T> data(get_col(col));
		for(MatrixColIter<T> i = data.begin(); i != data.end(); ++i) {
			vec.push_back(*i);
		}
	}

};

template<class T>
void Matrix<T>::resize(int num_rows, int num_cols, T def) {
	if(num_rows == 0 && num_cols == 0) {
		m_num_rows = num_rows;
		m_num_cols = num_cols;
		m_mat.resize(0);
		return;
	}
	if(num_cols != m_num_cols) {
		//ensure pointer validity after insert
		//cerr << "will rezize" << endl;
		vector<T> new_mat(num_cols * num_rows, def);
		vec_iter targ = new_mat.begin();
		for(vec_iter i = m_mat.begin();
		    i < m_mat.end();
		    i += m_num_cols) {
//			cerr << "copy at " << int(i-m_mat.begin()) << endl;
			copy(i, i + m_num_cols, targ);
			targ += num_cols;
		}
		m_num_cols = num_cols;
		m_mat = new_mat;
	} else {
		m_mat.resize(num_cols * num_rows, def);
	}
	m_num_rows = num_rows;
}

template<class T>
void print_csv_matrix(const Matrix<T> &mat, ostream &out, char delim = '\t')
{
	for(int r = 0; r < mat.row_size(); r++) {
		for(int c = 0; c < mat.col_size(); c++) {
			out << mat[r][c] << delim;
		}
		out << "\n";
	}
}
//matrix[i][j]

#endif // stdalg_ds_Matrix_h
