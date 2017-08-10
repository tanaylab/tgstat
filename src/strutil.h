#ifndef util_strutil_h
#define util_strutil_h 1

#include <vector>
#include <string>
#include "BufferedFile.h"

// If you know how many fields you're expecting, pass it via estimated_num_fields and split_line() will run twice faster

// split_line_by_space_chars reads and splits the line to fields separated by tab or space chars.
// Unlike other functions below split_line_by_space_chars eats up consequent space chars without creating empty fields.
// Returns the number of lines read (including the empty lines).
int split_line_by_space_chars(istream &in, vector<string> &fields, int estimated_num_fields = 1);
int split_line_by_space_chars(BufferedFile &in, vector<string> &fields, int estimated_num_fields = 1);
int split_line(istream &in, vector<string> &fields, char delim = '\t', int estimated_num_fields = 1);
int split_line(BufferedFile &in, vector<string> &fields, char delim, int estimated_num_fields = 1);
void split_line(const string &s, vector<string> &fields, char delim = '\t');
void split_line(const string &s, vector<float> &fields, char delim = '\t');
void split_line(const string &s, vector<int> &fields, char delim = '\t');
int get_one_field(istream &in, string &field, char delim, int num, bool eat_line);
void read_int_table(istream &in, int width, vector<vector<int> > &data);
void read_float_table(istream &in, int width, vector<vector<float> > &data);
void read_string_table(istream &in, int width, vector<vector<string> > &data);
int count_match(const string &targ, const string &mot);
void read_float_table_with_rowname(istream &in, vector<vector<float> > &data, vector<string> &row_name, int with_header, int subst_nas = 0, float na_value=0);

#endif // util_strutil_h
