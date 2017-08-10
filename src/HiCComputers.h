/*
 * Computer2D.h
 *
 *  Created on: May 21, 2012
 *      Author: eitany
 */

#ifndef HIC_COMPUTERS_H_
#define HIC_COMPUTERS_H_

#include <vector>

#include "Rectangle.h"
#include "BufferedFile.h"
#include "GenomeTrackSparse.h"
#include "Computer2D.h"

#include "Matrix.h"

using namespace std;

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// AreaComputer2D
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// returns const 1 which, when used with weighted sum results in the area
class AreaComputer2D: public Computer2D
{
public:
	AreaComputer2D(const char *trackdb_path, int chromid1, int chromid2) : Computer2D(CT2_AREA, trackdb_path, chromid1, chromid2) {}
    double compute(const Rectangle &rectangle) {
        return 1;
    }
    double compute(const Rectangle &rectangle, const DiagonalBand &band) {
    	if (!band.do_intersect(rectangle))
    		return 0;
    	if (band.do_contain(rectangle))
    		return 1;
    	Rectangle shrinked_rect(rectangle);
    	band.shrink2intersected(shrinked_rect);
    	return band.intersected_area(rectangle) / (double)rectangle.area();
    }
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// PotentialComputer2D
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class PotentialComputer2D: public Computer2D
{
public:
	PotentialComputer2D(const char *trackdb_path, int chromid1, int chromid2) : Computer2D(CT2_POTENTIAL, trackdb_path, chromid1, chromid2), m_loaded(false) {}
    double compute(const Rectangle &rectangle) { return compute(rectangle, NULL); }
    double compute(const Rectangle &rectangle, const DiagonalBand &band) { return compute(rectangle, &band); }

	void serialize(BufferedFile& bfile);
	void unserialize(BufferedFile& bfile);

	void dump();

    // filenames are relative to trackdb
    void set_fend_track(const char* track_fn1, const char* track_fn2);
protected:
    bool m_loaded;

    // need to serialize these:
    string m_track_fn1, m_track_fn2;

    // this is not serialized
    GenomeTrackSparse m_track1, m_track2;

    double compute(const Rectangle &rectangle, const DiagonalBand *band);
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// TechnicalComputer2D
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class TechnicalComputer2D: public Computer2D
{
public:
	TechnicalComputer2D(const char *trackdb_path, int chromid1, int chromid2)
        : Computer2D(CT2_TECHNICAL, trackdb_path, chromid1, chromid2), m_loaded(false), m_dim(0), m_prior(0), m_track1(0), m_track2(0) {}
   ~TechnicalComputer2D();
   double compute(const Rectangle &rectangle) { return compute(rectangle, NULL); }
   double compute(const Rectangle &rectangle, const DiagonalBand &band) { return compute(rectangle, &band); }

	void serialize(BufferedFile& bfile);
	void unserialize(BufferedFile& bfile);

	void dump();

    void set_prior(double prior);

    // filenames are relative to trackdb path
    void add_bias(const char* track_fn1, const char* track_fn2, const Matrix<double> matrix);
    void clear_biases();

protected:
    bool m_loaded;

    int m_dim;

    // prior trans contact
    double m_prior;

    // bias lookup tracks
    vector<string> m_track_fn1;
    vector<string> m_track_fn2;

    // bias matrices
    vector< Matrix<double> > m_matrix;

    // this is not serialized
    GenomeTrackSparse* m_track1;
    GenomeTrackSparse* m_track2;

    double compute(const Rectangle &rectangle, const DiagonalBand *band);
};


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// TestComputer2D
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Returns the Manhattan distance to the center of the intersected rectangle. This computer is used for testing.
class TestComputer2D: public Computer2D {
public:
	TestComputer2D(const char *trackdb_path, int chromid1, int chromid2) : Computer2D(CT2_TEST, trackdb_path, chromid1, chromid2) {}
    double compute(const Rectangle &rectangle) { return (rectangle.x1 + rectangle.x2 + rectangle.y1 + rectangle.y2) % 10000000; }
    double compute(const Rectangle &rectangle, const DiagonalBand &band) { return (rectangle.x1 + rectangle.x2 + rectangle.y1 + rectangle.y2 + band.d1 + band.d2) % 10000000; }
};

#endif
