#ifndef COMPUTER2D_H_
#define COMPUTER2D_H_

#include "BufferedFile.h"
#include "DiagonalBand.h"
#include "Rectangle.h"

using namespace std;

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Abstract computer class
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/* USAGE NOTES:
 1. When the class references 1D tracks it uses the trackdb_path as a basedir.
 2. Before calling unserialize you must call createComputer2D on the bfile to read the first byte which contains the type.
*/

class Computer2D {
public:
	typedef enum { CT2_AREA = 0, CT2_POTENTIAL, CT2_TECHNICAL, CT2_TEST  } Computer2DType;

	Computer2D(Computer2DType type, const char *trackdb_path, int chromid1, int chromid2)
        : m_type(type), m_trackdb_path(trackdb_path), m_chromid1(chromid1), m_chromid2(chromid2) {}
    virtual ~Computer2D() {}

    // get computer type
    Computer2DType get_type() { return m_type; }
    int get_chromid1() { return m_chromid1; }
    int get_chromid2() { return m_chromid2; }

    // used for debug, outputs object to screen
	virtual void dump() {}

    // compute value for rectangle
    virtual double compute(const Rectangle &rectangle) = 0;

    // compute value for rectangle
    virtual double compute(const Rectangle &rectangle, const DiagonalBand &band) = 0;

    // class functions
    static Computer2D* unserializeComputer2D(BufferedFile& bfile, const char *trackdb_path, int chromid1, int chromid2);
    static void serializeComputer2D(BufferedFile& bfile, Computer2D* computer);

protected:
    Computer2DType m_type;
    const char* m_trackdb_path;
    int m_chromid1, m_chromid2;

    // read/write computer
    // do not call this function directly, instead use serializeComputer2D and unserializeComputer2D
	virtual void serialize(BufferedFile& bfile) {}
	virtual void unserialize(BufferedFile& bfile) {}
};


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#endif
