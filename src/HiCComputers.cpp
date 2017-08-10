#include "TGLException.h"
#include "GenomeTrack2D.h"
#include "HiCComputers.h"

#include <fstream>

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Utility structure and functions
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

struct comp_interval_start : public binary_function<GInterval, GInterval, bool> {
    bool operator()(const GInterval x, const GInterval y) {
        if (x.chromid != y.chromid)
            return (x.chromid < y.chromid);
        else
            return (x.start < y.start);
    }
};

struct comp_interval_end : public binary_function<GInterval, GInterval, bool> {
    const bool operator()(const GInterval x, const GInterval y) {
        if (x.chromid != y.chromid)
            return (x.chromid < y.chromid);
        else
            return (x.end < y.end);
    }
};

static void binary_search(const GIntervals& ints, int chromid, int64_t coord, bool is_start, int& index)
{
    GIntervals::const_iterator it;
    GInterval interval(chromid, coord, coord, 1);
    if (is_start)
        it = upper_bound (ints.begin(), ints.end(), interval, comp_interval_end());
    else
        it = lower_bound (ints.begin(), ints.end(), interval, comp_interval_start());
    index = it - ints.begin();
}

static void write_string(BufferedFile& bfile, string str)
{
	size_t size = str.size();
    bfile.write(&size, sizeof(size));
    if (bfile.write(&*str.begin(), sizeof(char)*size) != size)
		TGLError("Writing string failed, file: %s", bfile.file_name().c_str());
}

static void read_string(BufferedFile& bfile, string& str)
{
	size_t size;
    bfile.read(&size, sizeof(size));

    str.resize(size);
    if (bfile.read(&*str.begin(), sizeof(char)*size) != size)
		TGLError("Reading string failed, file: %s", bfile.file_name().c_str());
}

static void write_matrix(BufferedFile& bfile, Matrix<double>& matrix)
{
    int num_rows = matrix.row_size();
    int num_cols = matrix.col_size();
    bfile.write(&num_rows, sizeof(num_rows));
    bfile.write(&num_cols, sizeof(num_cols));

    vector<double>& vec = matrix.get_vector();
    if (bfile.write(&*vec.begin(), sizeof(double)*vec.size()) != sizeof(double)*vec.size())
		TGLError("Writing matrix failed, file: %s", bfile.file_name().c_str());
}

static void read_matrix(BufferedFile& bfile, Matrix<double>& matrix)
{
    int num_rows;
    int num_cols;
    bfile.read(&num_rows, sizeof(num_rows));
    bfile.read(&num_cols, sizeof(num_cols));
    matrix.resize(num_rows, num_cols);

    vector<double>& vec = matrix.get_vector();
    if (bfile.read(&*vec.begin(), sizeof(double)*vec.size()) != sizeof(double)*vec.size())
		TGLError("Reading matrix failed, file: %s", bfile.file_name().c_str());
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// PotentialComputer2D
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double PotentialComputer2D::compute(const Rectangle &rectangle, const DiagonalBand *band)
{
    if (!m_loaded)
    {
        string track_fn1 = string(m_trackdb_path) + "/" + m_track_fn1;
        string track_fn2 = string(m_trackdb_path) + "/" + m_track_fn2;
        m_track1.init_read(track_fn1.c_str(), m_chromid1);
        m_track2.init_read(track_fn2.c_str(), m_chromid2);
        m_loaded = true;
    }

    GInterval interval1(m_chromid1, rectangle.x1, rectangle.x2, 1);
    GInterval interval2(m_chromid2, rectangle.y1, rectangle.y2, 1);

    // cout << "x1=" << rectangle.x1 << " y1=" << rectangle.y1 << " x2=" << rectangle.x2 << " y2=" << rectangle.y2 << endl;

    int start1, end1, start2, end2;

    const GIntervals& intervals1 = m_track1.get_intervals();
    const GIntervals& intervals2 = m_track2.get_intervals();
    if ((intervals1.size() == 0) || (intervals2.size() == 0))
        return 0;

    // find indices matching the query rectantgle
    binary_search(intervals1, m_chromid1, interval1.start, true, start1);
    binary_search(intervals1, m_chromid1, interval1.end, false, end1);
    binary_search(intervals2, m_chromid2, interval2.start, true, start2);
    binary_search(intervals2, m_chromid2, interval2.end, false, end2);

    // cout << "indices: start1=" << start1 << ", end1=" << end1 << ", start2=" << start2 << ", end2=" << end2 << endl;

//     cout << "start1: " << intervals1[start1].start << " - " << intervals1[start1].end << endl;
//     cout << "end1: " << intervals1[end1-1].start << " - " << intervals1[end1-1].end << endl;
//     cout << "start2: " << intervals2[start2].start << " - " << intervals2[start2].end << endl;
//     cout << "end2: " << intervals2[end2-1].start << " - " << intervals2[end2-1].end << endl;

    double result = 0;
    for (int i1 = start1; i1 < end1; i1++) {
    for (int i2 = start2; i2 < end2; i2++) {
    	if (band && band->do_contain(Rectangle(intervals1[i1].start, intervals2[i2].start, intervals1[i1].end, intervals2[i2].end)))
            continue;
        result++;
    } }

    // divide by area
    double area = (rectangle.x2 - rectangle.x1) * (rectangle.y2 - rectangle.y1);
    result = result / area;

    return result;
}

void PotentialComputer2D::serialize(BufferedFile& bfile)
{
    write_string(bfile, m_track_fn1);
    write_string(bfile, m_track_fn2);
}

void PotentialComputer2D::unserialize(BufferedFile& bfile)
{
    read_string(bfile, m_track_fn1);
    read_string(bfile, m_track_fn2);
}

void PotentialComputer2D::dump()
{
    cout << "m_type: " << m_type << endl;
    cout << "m_chromid1: " << m_chromid1 << endl;
    cout << "m_chromid2: " << m_chromid2 << endl;
    cout << "m_track_fn1: " << m_track_fn1 << endl;
    cout << "m_track_fn2: " << m_track_fn2 << endl;
}

void PotentialComputer2D::set_fend_track(const char* track_fn1, const char* track_fn2)
{
    m_track_fn1 = string(track_fn1);
    m_track_fn2 = string(track_fn2);
    m_loaded = false;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// TechnicalComputer2D
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


TechnicalComputer2D::~TechnicalComputer2D()
{
    delete[] m_track1;
    delete[] m_track2;
}

static bool file_exists(const char * filename)
{
    std::fstream f;
    f.open(filename);
    return (f.is_open());
}

double TechnicalComputer2D::compute(const Rectangle &rectangle, const DiagonalBand *band)
{
//	printf("tech pcomputer called\n");

	if ( (rectangle.x1 >= rectangle.x2) || (rectangle.y1 >= rectangle.y2) )
        return 0;

    if (m_dim == 0)
        TGLError("Assertion: dim must be >0");

    if ( ((int)m_track_fn1.size() != (int)m_dim) || ((int)m_track_fn2.size() != m_dim) || ((int)m_matrix.size() != m_dim))
        TGLError("Assertion: number of correction tracks (d1=%d, d2=%d) and number of correction matrices (%d) must be equal to %d",
                 m_track_fn1.size(), m_track_fn2.size(), m_matrix.size(), m_dim);

    if (!m_loaded)
    {
        delete[] m_track1;
        delete[] m_track2;
        m_track1 = new GenomeTrackSparse[m_dim];
        m_track2 = new GenomeTrackSparse[m_dim];

        for (int i = 0; i < m_dim; i++) {
            // track1
            string track_fn1 = string(m_trackdb_path) + "/" + m_track_fn1[i];
            if ( !file_exists(track_fn1.c_str()))
                TGLError("File does not exist: %s", track_fn1.c_str());
            m_track1[i].init_read(track_fn1.c_str(), m_chromid1);
            // track2
            string track_fn2 = string(m_trackdb_path) + "/" + m_track_fn2[i];
            if ( !file_exists(track_fn2.c_str()))
                TGLError("File does not exist: %s", track_fn2.c_str());
            m_track2[i].init_read(track_fn2.c_str(), m_chromid2);
        }
        m_loaded = true;
    }

    GInterval interval1(m_chromid1, rectangle.x1, rectangle.x2, 1);
    GInterval interval2(m_chromid2, rectangle.y1, rectangle.y2, 1);

    int start1=0, end1=0, start2=0, end2=0;

    const GIntervals& gintervals1 = m_track1[0].get_intervals();
    const GIntervals& gintervals2 = m_track2[0].get_intervals();

    for (int i = 0; i < m_dim; i++) {
        const GIntervals& intervals1 = m_track1[i].get_intervals();
        const GIntervals& intervals2 = m_track2[i].get_intervals();

        int tstart1, tend1, tstart2, tend2;

        // find indices matching the query rectantgle
        binary_search(intervals1, m_chromid1, interval1.start, true, tstart1);
        binary_search(intervals1, m_chromid1, interval1.end, false, tend1);
        binary_search(intervals2, m_chromid2, interval2.start, true, tstart2);
        binary_search(intervals2, m_chromid2, interval2.end, false, tend2);

        if (i == 0)
        {
            start1 = tstart1;
            start2 = tstart2;
            end1 = tend1;
            end2 = tend2;
        }

        if ( (start1 != tstart1) || (start2 != tstart2) || (end1 != tend1) || (end2 != tend2) )
            TGLError("model tracks must have identical structure");
    }
    if ((start1 == end1) || (start2 == end2))
        return 0;

    // cout << "start1=" << start1 << " end1=" << end1 << " start2=" << start2 << " end2=" << end2 << endl;
    //printf("prior=%e\n", m_prior);

    double result = 0;
    for (int i1 = start1; i1 < end1; i1++) {
    for (int i2 = start2; i2 < end2; i2++) {
    	if (band && !band->do_contain(Rectangle(gintervals1[i1].start, gintervals2[i2].start, gintervals1[i1].end, gintervals2[i2].end)))
            continue;
        double value = m_prior;
        for(unsigned int mi = 0; mi < m_matrix.size(); mi++) {
            int f1 = (int)m_track1[mi].get_vals()[i1] - 1;
            int f2 = (int)m_track2[mi].get_vals()[i2] - 1;
            value *= m_matrix[mi][f1][f2];
            // printf("f1=%d, f2=%d, factor=%.10f   ////", f1, f2, m_matrix[mi][f1][f2]);
        }
        // cout << endl;
        result += value;
    } }

    double area = (rectangle.x2 - rectangle.x1) * (rectangle.y2 - rectangle.y1);
    result /= area;

    return result;
}

void TechnicalComputer2D::serialize(BufferedFile& bfile)
{
    // write dimension
	bfile.write(&m_dim, sizeof(m_dim));

    // write prior
	bfile.write(&m_prior, sizeof(m_prior));

    for (int i = 0; i < m_dim; i++)
    {
        write_string(bfile, m_track_fn1[i]);
        write_string(bfile, m_track_fn2[i]);
        write_matrix(bfile, m_matrix[i]);
    }
}

void TechnicalComputer2D::unserialize(BufferedFile& bfile)
{
    // read dimension
	bfile.read(&m_dim, sizeof(m_dim));

    // write prior
	bfile.read(&m_prior, sizeof(m_prior));

    // resize structs
    m_track_fn1.resize(m_dim);
    m_track_fn2.resize(m_dim);
    m_matrix.resize(m_dim);

    for (int i = 0; i < m_dim; i++)
    {
        read_string(bfile, m_track_fn1[i]);
        read_string(bfile, m_track_fn2[i]);
        read_matrix(bfile, m_matrix[i]);
    }
}

void TechnicalComputer2D::dump()
{
    cout << "m_type: " << m_type << endl;
    cout << "m_chromid1: " << m_chromid1 << endl;
    cout << "m_chromid2: " << m_chromid2 << endl;

    cout << "m_dim: " << m_dim << endl;
    for (unsigned int i = 0; i < m_track_fn1.size(); i++) {
        cout << "m_track_fn1: " << m_track_fn1[i] << endl;
        cout << "m_track_fn2: " << m_track_fn2[i] << endl;
        cout << "m_matrix: num_cols=" << m_matrix[i].col_size() << " num_rows=" << m_matrix[i].row_size() << endl;
//         for (unsigned int row = 0; row < m_matrix[i].row_size(); row++) {
//             for (unsigned int col = 0; col < m_matrix[i].col_size(); col++) {
//                 cout << m_matrix[i][col][row] << " ";
//             }
//             cout << endl;
//         }
    }
}

void TechnicalComputer2D::set_prior(double prior)
{
    m_prior = prior;
    m_loaded = false;
}

void TechnicalComputer2D::add_bias(const char* track_fn1, const char* track_fn2, const Matrix<double> matrix)
{
    m_track_fn1.push_back(track_fn1);
    m_track_fn2.push_back(track_fn2);
    m_matrix.push_back(matrix);

    m_dim++;

    m_loaded = false;
}

void TechnicalComputer2D::clear_biases()
{
    m_track_fn1.resize(0);
    m_track_fn2.resize(0);
    m_matrix.resize(0);

    m_dim = 0;

    m_loaded = false;
}
