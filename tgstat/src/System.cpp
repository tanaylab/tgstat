#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "BufferedFile.h"
#include "System.h"
#include "port.h"
#include "strutil.h"

size_t get_unique_mem_usage(pid_t pid)
{
	size_t mem_usage = 0;
	BufferedFile bf;
	char filename[100];
	vector<string> fields;

	sprintf(filename, "/proc/%ld/smaps", (long)pid);

	// count only private dirty bytes under heap section
	if (!bf.open(filename, "r")) {
		bool is_heap = false;

		while (1) {
			split_line_by_space_chars(bf, fields, 2);

			if (fields.empty()) 
				break;

			if (fields.size() == 6) {
				if (is_heap) {
					if (!fields[5].empty() && fields[5] != "[heap]")
						break;
				} else if (fields[5] == "[heap]") 
					is_heap = true;
			} else if (is_heap && fields.size() == 3 && fields[0] == "Private_Dirty:") {
				char *endptr;
				long num;

				num = strtol(fields[1].c_str(), &endptr, 10);
				if (!*endptr)
					mem_usage += num;
			}
		}
	}
	return mem_usage * 1024;
}

